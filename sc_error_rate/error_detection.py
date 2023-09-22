from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pysam
from Bio import SeqIO, Seq
from dataclasses_json import DataClassJsonMixin
from tqdm import tqdm

from sc_error_rate.paths import split_by_cells_folder_path, reference_genome_folder_path, \
    detected_errors_data_folder_path

DNA_ALPHABET = ('A', 'C', 'G', 'T')


@dataclass
class CountedSymbol(DataClassJsonMixin):
    symbol: str
    count: int


@dataclass
class DetectedError(DataClassJsonMixin):
    chromosome_name: str
    position: int
    reference_nucleotide: str
    detected_nucleotide: str
    num_unique_umi: int
    read_counts: list[CountedSymbol]
    consensus_read_counts: list[CountedSymbol]


@dataclass
class CellErrorStatistics(DataClassJsonMixin):
    barcode: str

    num_positions_with_sufficient_coverage: int = 0
    num_positions_with_dna_mutations_detected: int = 0
    num_positions_with_transcription_error_detected: int = 0

    num_consensus_reads_total: int = 0
    num_consensus_reads_with_dna_mutation: int = 0
    num_consensus_reads_with_transcription_error: int = 0

    transcription_error_rate_per_position: Optional[float] = None
    transcription_error_rate_per_consensus_read: Optional[float] = None
    mutation_error_rate_per_position: Optional[float] = None
    mutation_error_rate_per_consensus_read: Optional[float] = None

    transcription_errors: list[DetectedError] = field(default_factory=list)
    mutation_errors: list[DetectedError] = field(default_factory=list)

    def compute_derived_statistics(self):
        if self.num_positions_with_sufficient_coverage > 0:
            self.transcription_error_rate_per_position = self.num_positions_with_transcription_error_detected / self.num_positions_with_sufficient_coverage
            self.transcription_error_rate_per_consensus_read = self.num_positions_with_transcription_error_detected / self.num_consensus_reads_total

            self.mutation_error_rate_per_position = self.num_positions_with_dna_mutations_detected / self.num_positions_with_sufficient_coverage
            self.mutation_error_rate_per_consensus_read = self.num_consensus_reads_with_dna_mutation / self.num_consensus_reads_total


@dataclass
class AlignedReadWithPosition:
    read: pysam.AlignedSegment
    position: int


def get_consensus_read(aligned_reads_with_position: list[AlignedReadWithPosition]) -> Optional[str]:
    """
    Performs majority voting for a nucleotide at given position.
    :param aligned_reads_with_position: list of AlignedReadWithPosition, each storing the pysam.AlignedSegment object
    and position in the segment to which the reference sequence is aligned.
    :return: Optional symbol from the DNA_ALPHABET. If there are 3 and more different nucleotides detected, or the
     voting result in a tie, None is returned.
    """
    symbol_occurrences = count_symbol_occurrences(
        [aligned_read.read.query_sequence[aligned_read.position] for aligned_read in aligned_reads_with_position])
    if symbol_occurrences[2].count > 0:
        return None  # more than 2 different bases at single position
    return symbol_occurrences[0].symbol if symbol_occurrences[0].count > len(aligned_reads_with_position) / 2 else None


def count_symbol_occurrences(nucleotides: list[str]) -> list[CountedSymbol]:
    """
    Counts the occurrences of nucleotides in the input list.
    :param nucleotides: List of nucleotides, which are assumed to be upper-case symbols from the DNA_ALPHABET ('A', 'C',
     'G' or 'T') .
    :return: List of CountedSymbol object. Each CountedSymbol represent one symbol from DNA_ALPHABET and the count of
     its occurrences in the nucleotides list. The output list is sorted by the symbol count in descending manner
     (i.e., first element in the list corresponds to the symbol with the most occurrences).
    """
    nucleotide_counts = {symbol: 0 for symbol in DNA_ALPHABET}
    for nucleotide in nucleotides:
        nucleotide_counts[nucleotide] += 1
    counted_symbols_list = [CountedSymbol(symbol=symbol, count=count) for symbol, count in nucleotide_counts.items()]
    return sorted(counted_symbols_list, key=lambda x: x.count, reverse=True)


# %%
def compute_cell_error_statistics(bam_file_path: Path,
                                  reference_genome: dict[str, Seq],
                                  min_reads_per_umi=5,
                                  min_consensus_reads=5) -> CellErrorStatistics:
    """
    Detects transcriptional errors and mutations in the reads from single cell.
    :param bam_file_path: Path the a .bam file, which contains reads from a single cell.
    :param reference_genome: Dictionary of chromosome names and corresponding sequences.
    :param min_reads_per_umi: Minimum amount of reads with the same UMI needed to create a 'consensus read'.
    :param min_consensus_reads: Minimum amount of consensus reads needed to detect an error.
    :return: CellErrorStatistics object, which store the information about detected errors and additional statistics
    about coverage of the reference sequence.
    """
    samfile = pysam.AlignmentFile(bam_file_path, "rb")
    cell_error_statistics = CellErrorStatistics(barcode=bam_file_path.stem)

    for chromosome_number in list(range(1, 20)) + ['X', 'Y']:  # skipping MT chromosome for now
        chromosome_name = f"chr{chromosome_number}"
        for pileup_column in samfile.pileup(chromosome_name):
            if pileup_column.get_num_aligned() < min_reads_per_umi * min_consensus_reads:
                continue

            mapping_qualities = pileup_column.get_mapping_qualities()
            aligned_reads = [pileup_read.alignment for pileup_read, mapping_quality in
                             zip(pileup_column.pileups, mapping_qualities) if mapping_quality >= 255]
            query_positions = [position for position, mapping_quality in
                               zip(pileup_column.get_query_positions(), mapping_qualities) if mapping_quality >= 255]

            aligned_reads_by_umi: dict[str, list[AlignedReadWithPosition]] = defaultdict(list)
            for read, position in zip(aligned_reads, query_positions):
                aligned_reads_by_umi[read.get_tag('UB')].append(AlignedReadWithPosition(read=read, position=position))

            if len(aligned_reads_by_umi) < min_consensus_reads:
                continue

            consensus_reads: list[str] = []
            for aligned_reads_with_position in aligned_reads_by_umi.values():
                if len(aligned_reads_with_position) < min_reads_per_umi:
                    continue
                consensus_read = get_consensus_read(aligned_reads_with_position)
                if consensus_read is not None:
                    consensus_reads.append(consensus_read)

            if len(consensus_reads) < min_consensus_reads:
                continue

            consensus_read_counts = count_symbol_occurrences(consensus_reads)
            num_consensus_reads = len(consensus_reads)
            reference_symbol = reference_genome[chromosome_name][pileup_column.reference_pos]

            if num_consensus_reads < min_consensus_reads:
                continue  # not enough consensus reads

            cell_error_statistics.num_positions_with_sufficient_coverage += 1
            cell_error_statistics.num_consensus_reads_total += num_consensus_reads

            if ((consensus_read_counts[0].symbol == reference_symbol) and
                    (consensus_read_counts[0].count == num_consensus_reads)):
                continue  # all consensus reads agree with the reference

            detected_error = DetectedError(chromosome_name=chromosome_name,
                                           position=pileup_column.reference_pos,
                                           reference_nucleotide=reference_symbol,
                                           detected_nucleotide=consensus_read_counts[0].symbol if consensus_read_counts[
                                                                                                      0].symbol != reference_symbol else
                                           consensus_read_counts[1].symbol,
                                           num_unique_umi=len(aligned_reads_by_umi),
                                           read_counts=count_symbol_occurrences(
                                               [read.query_sequence[position] for read, position in
                                                zip(aligned_reads, query_positions)]),
                                           consensus_read_counts=consensus_read_counts)
            if ((consensus_read_counts[0].symbol == reference_symbol) and
                    (consensus_read_counts[0].count == num_consensus_reads - 1)):
                # exactly one consensus read differs from the reference sequence,
                # which we classify as transcriptional error
                cell_error_statistics.transcription_errors.append(detected_error)
                cell_error_statistics.num_positions_with_transcription_error_detected += 1
                cell_error_statistics.num_consensus_reads_with_transcription_error += 1
            else:
                # more than one consensus read differs from the reference - we classify this as DNA mutation / damage
                cell_error_statistics.mutation_errors.append(detected_error)
                cell_error_statistics.num_positions_with_dna_mutations_detected += 1
                cell_error_statistics.num_consensus_reads_with_dna_mutation += sum(
                    [x.count for x in consensus_read_counts if x.symbol != reference_symbol])
    cell_error_statistics.compute_derived_statistics()
    return cell_error_statistics


def process_cell_files(input_folder: Path,
                       output_folder: Path,
                       reference_genome_fasta_file_path: Path) -> None:
    """
    Iterates through all .bam file in given folder. For each file, computes CellErrorStatistics and write them to
    JSON file in the output folder (the name of the output file is the same as of the input file).
    :param input_folder: Input folder with .bam files, corresponding to individual cells.
    :param output_folder: Folder to which JSON files with results will be written.
    :param reference_genome_fasta_file_path: Path to fasta file with reference genome.
    :return: None
    """
    reference_genome_all_sequences = list(SeqIO.parse(open(reference_genome_fasta_file_path), 'fasta'))
    reference_genome: dict[str, Seq] = {record.name: record.seq for record in reference_genome_all_sequences
                                        if record.name.startswith('chr')}
    output_folder.mkdir(parents=True, exist_ok=True)
    bam_file_paths = [file_path for file_path in input_folder.iterdir() if file_path.suffix == '.bam']
    for file_path in tqdm(bam_file_paths, desc='Detecting sequence errors in the .bam files of individual cells'):
        cell_error_statistics = compute_cell_error_statistics(file_path, reference_genome)
        with open(output_folder / f'{cell_error_statistics.barcode}.json', 'w') as output_file:
            output_file.write(cell_error_statistics.to_json())


if __name__ == "__main__":
    process_cell_files(input_folder=split_by_cells_folder_path / 'DR1_old',
                       output_folder=detected_errors_data_folder_path / 'DR1_old_255_quality',
                       reference_genome_fasta_file_path=reference_genome_folder_path / 'genome.fa')

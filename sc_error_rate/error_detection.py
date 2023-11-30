import argparse
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from threading import Thread
from typing import Optional

import pysam

from Bio import SeqIO, Seq
from dataclasses_json import DataClassJsonMixin
from tqdm import tqdm

DNA_ALPHABET = ('A', 'C', 'G', 'T')


@dataclass
class CountedSymbol(DataClassJsonMixin):
    symbol: str
    count: int


@dataclass
class MutationError(DataClassJsonMixin):
    chromosome_name: str
    position: int
    reference_nucleotide: str
    detected_nucleotide: str
    num_unique_umi: int
    read_counts: list[CountedSymbol]
    consensus_read_counts: list[CountedSymbol]


@dataclass
class TranscriptionError(DataClassJsonMixin):
    chromosome_name: str
    position: int
    reference_nucleotide: str
    majority_nucleotide: str
    error_nucleotide: str
    error_umi: str
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

    transcription_errors: list[TranscriptionError] = field(default_factory=list)
    mutation_errors: list[MutationError] = field(default_factory=list)

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
    return symbol_occurrences[0].symbol if symbol_occurrences[0].count == len(aligned_reads_with_position) else None
    # return symbol_occurrences[0].symbol if symbol_occurrences[0].count > len(aligned_reads_with_position) / 2 else None


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
                                  min_consensus_reads=5,
                                  record_mutations=False) -> CellErrorStatistics:
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
    chromosomes_of_interest = [reference_sequence['SN'] for reference_sequence in samfile.header.to_dict()['SQ'] if
                               reference_sequence['SN'].startswith('chr') and reference_sequence['SN'] != 'chrM']
    assert chromosomes_of_interest

    for chromosome_name in chromosomes_of_interest:  # skipping MT chromosome for now
        for pileup_column in samfile.pileup(chromosome_name):
            if pileup_column.get_num_aligned() < min_reads_per_umi * min_consensus_reads:
                continue

            # mapping_qualities = pileup_column.get_mapping_qualities()
            # aligned_reads = [pileup_read.alignment for pileup_read, mapping_quality in
            #                  zip(pileup_column.pileups, mapping_qualities) if mapping_quality >= 255]
            # query_positions = [position for position, mapping_quality in
            #                    zip(pileup_column.get_query_positions(), mapping_qualities) if mapping_quality >= 255]
            #
            # aligned_reads_by_umi: dict[str, list[AlignedReadWithPosition]] = defaultdict(list)
            # for read, position in zip(aligned_reads, query_positions):
            #     aligned_reads_by_umi[read.get_tag('UB')].append(AlignedReadWithPosition(read=read, position=position))
            aligned_reads_by_umi: dict[str, list[AlignedReadWithPosition]] = defaultdict(list)
            for pileup_read in pileup_column.pileups:
                if pileup_read.query_position is not None:
                    read = pileup_read.alignment
                    if read.mapping_quality >= 255:
                        aligned_reads_by_umi[read.get_tag('UB')].append(
                            AlignedReadWithPosition(read=read, position=pileup_read.query_position))

            if len(aligned_reads_by_umi) < min_consensus_reads:
                continue

            consensus_reads: list[str] = []
            umis_by_symbol = {symbol: [] for symbol in DNA_ALPHABET}
            for umi, aligned_reads_with_position in aligned_reads_by_umi.items():
                if len(aligned_reads_with_position) < min_reads_per_umi:
                    continue
                consensus_read = get_consensus_read(aligned_reads_with_position)
                if consensus_read is not None:
                    consensus_reads.append(consensus_read)
                    umis_by_symbol[consensus_read].append(umi)

            if len(consensus_reads) < min_consensus_reads:
                continue

            num_consensus_reads = len(consensus_reads)
            if num_consensus_reads < min_consensus_reads:
                continue  # not enough consensus reads

            consensus_read_counts = count_symbol_occurrences(consensus_reads)
            reference_symbol = reference_genome[chromosome_name][pileup_column.reference_pos]

            cell_error_statistics.num_positions_with_sufficient_coverage += 1
            cell_error_statistics.num_consensus_reads_total += num_consensus_reads

            if consensus_read_counts[0].count == num_consensus_reads and consensus_read_counts[
                0].symbol == reference_symbol:
                # all UMIs agree with the reference genome
                pass
            if consensus_read_counts[0].count == num_consensus_reads - 1:
                assert consensus_read_counts[1].count == 1
                assert len(umis_by_symbol[consensus_read_counts[1].symbol]) == 1
                # exactly one UMI differs from others - we consider this to be a transcription error
                cell_error_statistics.num_positions_with_transcription_error_detected += 1
                cell_error_statistics.num_consensus_reads_with_transcription_error += 1
                cell_error_statistics.transcription_errors.append(
                    TranscriptionError(chromosome_name=chromosome_name,
                                       position=pileup_column.reference_pos,
                                       reference_nucleotide=reference_symbol,
                                       majority_nucleotide=consensus_read_counts[0].symbol,
                                       error_nucleotide=consensus_read_counts[1].symbol,
                                       error_umi=umis_by_symbol[consensus_read_counts[1].symbol][0],
                                       num_unique_umi=len(aligned_reads_by_umi),
                                       read_counts=count_symbol_occurrences(
                                           [read.query_sequence[position] for read, position in
                                            zip(aligned_reads, query_positions)]),
                                       consensus_read_counts=consensus_read_counts
                                       )
                )

            elif consensus_read_counts[0].count != num_consensus_reads or consensus_read_counts[
                0].symbol != reference_symbol:
                # mutation error
                cell_error_statistics.num_positions_with_dna_mutations_detected += 1
                cell_error_statistics.num_consensus_reads_with_dna_mutation += sum(
                    [x.count for x in consensus_read_counts if x.symbol != reference_symbol])
                if record_mutations:
                    raise NotImplementedError  # TODO
            else:
                # all UMIs should agree with the reference genome
                assert consensus_read_counts[0].count == num_consensus_reads and consensus_read_counts[
                    0].symbol == reference_symbol

    cell_error_statistics.compute_derived_statistics()
    return cell_error_statistics


def process_cell_files(input_folder: Path,
                       output_folder: Path,
                       reference_genome_fasta_file_path: Path,
                       record_mutations=False) -> None:
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
        cell_error_statistics = compute_cell_error_statistics(file_path,
                                                              reference_genome,
                                                              record_mutations=record_mutations)
        with open(output_folder / f'{cell_error_statistics.barcode}.json', 'w') as output_file:
            output_file.write(cell_error_statistics.to_json())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_folder',
                        default='/cellfile/datapublic/jkoubele/leibniz_institute_data/split_by_cells/20201014_582_KLR/GEX/O_DR')
    parser.add_argument('--output_folder',
                        default='/cellfile/datapublic/jkoubele/transcription_error_analysis/computed_data/transcription_errors/FLI_O_DR')
    parser.add_argument('--reference_genome_fasta_file',
                        default='/cellfile/datapublic/jkoubele/STAR_2.7.11a/reference_genomes/GRCm38/GRCm38.primary_assembly.genome.fa')
    parser.add_argument('--record_mutations', action='store_true')
    args = parser.parse_args()
    print(f"{args=}")
    process_cell_files(input_folder=Path(args.input_folder),
                       output_folder=Path(args.output_folder),
                       reference_genome_fasta_file_path=Path(args.reference_genome_fasta_file),
                       record_mutations=args.record_mutations)

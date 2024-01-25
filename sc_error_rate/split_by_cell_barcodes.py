import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'  # solves weird error when importing numpy (and consequently e.g. pandas, biopython etc.) on cluster

from pathlib import Path
import argparse
import pandas as pd
import pysam
from tqdm import tqdm


def split_bam_file_by_cell_barcodes(bam_file_path: Path,
                                    cell_barcodes: set[str],
                                    output_folder: Path) -> None:
    """
    Split .bam file to multiple smaller ones, each containing reads from a single cell.
    :param bam_file_path: Path to a .bam file with reads.
    :param cell_barcodes: Set of cell barcodes. Reads with barcodes not present in cell_barcodes will be ignored.
    :param output_folder: Folder to which the resulting files will be written.
    :return: None.
    """
    output_folder.mkdir(parents=True, exist_ok=True)
    samfile_input = pysam.AlignmentFile(bam_file_path, "rb")
    output_file_names_by_barcode = {barcode: output_folder / f"{barcode}.bam" for barcode in cell_barcodes}
    output_samfiles_by_barcode = {barcode: pysam.AlignmentFile(file_path, "wb",
                                                               template=samfile_input)
                                  for barcode, file_path in output_file_names_by_barcode.items()}
    for read in tqdm(samfile_input, desc='Processing reads from the input .bam file'):
        if not (read.has_tag('CB') and read.has_tag('UB')):
            continue

        read_cell_barcode = read.get_tag('CB')
        umi = read.get_tag('UB')
        if len(read_cell_barcode) <= 5 or len(umi) <= 5:
            continue  # instead of CB or UB tag being None, it can be also '-', ' - ' or similar
        if read_cell_barcode in cell_barcodes:
            output_samfiles_by_barcode[read_cell_barcode].write(read)

    for samfile_output in output_samfiles_by_barcode.values():
        samfile_output.close()
    for file_name in tqdm(output_file_names_by_barcode.values(), desc='Indexing .bam files'):
        pysam.index(str(file_name))


def split_star_solo_output(star_solo_output_folder: Path, output_folder: Path) -> None:
    barcodes_df = pd.read_csv(
        star_solo_output_folder / Path('Solo.out/GeneFull/filtered/barcodes.tsv'),
        delimiter='\t',
        names=['barcode'])
    split_bam_file_by_cell_barcodes(
        bam_file_path=star_solo_output_folder / 'Aligned.sortedByCoord.out.bam',
        cell_barcodes=set(barcodes_df['barcode']),
        output_folder=output_folder)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam_input_folder',
                        default='/cellfile/datapublic/jkoubele/perturb_seq_cf/BAM/SRR21603868',
                        help='Folder containing .bam file to be processed.')
    parser.add_argument('--output_folder',
                        default='/cellfile/datapublic/jkoubele/perturb_seq_cf/split_by_cells/SRR21603868',
                        help='Output .bam files will be written here.')
    args = parser.parse_args()
    split_star_solo_output(star_solo_output_folder=Path(args.bam_input_folder),
                           output_folder=Path(args.output_folder))

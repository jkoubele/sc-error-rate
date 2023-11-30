from pathlib import Path

import pandas as pd
import pysam
from tqdm import tqdm

from sc_error_rate.paths import split_by_cells_folder_path, raw_data_folder_path


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

        read_barcode = read.get_tag('CB')
        if read_barcode in cell_barcodes:
            output_samfiles_by_barcode[read_barcode].write(read)

    for samfile_output in output_samfiles_by_barcode.values():
        samfile_output.close()
    for file_name in tqdm(output_file_names_by_barcode.values(), desc='Indexing .bam files'):
        pysam.index(str(file_name))


def split_aging_mouse_data_by_cell_barcodes(aging_mouse_data_folder_path: Path) -> None:
    barcodes_df = pd.read_csv(
        aging_mouse_data_folder_path / Path('outs/filtered_gene_bc_matrices/gencode.vM19/barcodes.tsv'),
        delimiter='\t',
        names=['barcode'])
    split_bam_file_by_cell_barcodes(bam_file_path=aging_mouse_data_folder_path / 'outs/possorted_genome_bam.bam',
                                    cell_barcodes=set(barcodes_df['barcode']),
                                    output_folder=split_by_cells_folder_path / aging_mouse_data_folder_path.name)


def split_dietary_restriction_mouse_data_by_cell_barcodes(dietary_restriction_mouse_data_folder_path: Path) -> None:
    barcodes_df = pd.read_csv(
        dietary_restriction_mouse_data_folder_path / Path('Solo.out/GeneFull/filtered/barcodes.tsv'),
        delimiter='\t',
        names=['barcode'])
    split_bam_file_by_cell_barcodes(
        bam_file_path=dietary_restriction_mouse_data_folder_path / 'Aligned.sortedByCoord.out.bam',
        cell_barcodes=set(barcodes_df['barcode']),
        output_folder=split_by_cells_folder_path / dietary_restriction_mouse_data_folder_path.name)


if __name__ == "__main__":
    # cellfile_data_path = Path('/cellfile/datapublic/jkoubele/sc_error_rate_data')
    input_folder_path = Path('/cellfile/datapublic/acherdi1/rnaspeed/datasets_mm/age_DR/res/young_AL5')
    # pysam.index(str(input_folder_path / 'Aligned.sortedByCoord.out.bam'))
    barcodes_df = pd.read_csv(
        input_folder_path / Path('Solo.out/GeneFull/filtered/barcodes.tsv'),
        delimiter='\t',
        names=['barcode'])
    split_bam_file_by_cell_barcodes(
        bam_file_path=input_folder_path / 'Aligned.sortedByCoord.out.bam',
        cell_barcodes=set(barcodes_df['barcode']),
        output_folder=Path('/cellfile/datapublic/jkoubele/DR_dataset/young_AL5'))

    input_folder_path_2 = Path('/cellfile/datapublic/acherdi1/rnaspeed/datasets_mm/age_DR/res/AL4_old')
    # pysam.index(str(input_folder_path / 'Aligned.sortedByCoord.out.bam'))
    barcodes_df = pd.read_csv(
        input_folder_path_2 / Path('Solo.out/GeneFull/filtered/barcodes.tsv'),
        delimiter='\t',
        names=['barcode'])
    split_bam_file_by_cell_barcodes(
        bam_file_path=input_folder_path_2 / 'Aligned.sortedByCoord.out.bam',
        cell_barcodes=set(barcodes_df['barcode']),
        output_folder=Path('/cellfile/datapublic/jkoubele/DR_dataset/AL4_old'))

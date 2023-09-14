from pathlib import Path
from tqdm import tqdm
import pysam
from sc_error_rate.paths import split_by_cells_folder_path, raw_data_folder_path
import pandas as pd


def split_bam_file_by_cell_barcodes(bam_file_path: Path,
                                    cell_barcodes: set[str],
                                    output_folder: Path) -> None:
    """
    :param bam_file_path: Path to a .bam file with reads.
    :param cell_barcodes: Set of cell barcodes. Reads with barcodes not present in cell_barcodes will be ignored.
    :param output_folder: Folder to which the resulting files will be written.
    :return: None.
    """
    output_folder.mkdir(parents=True, exist_ok=True)
    samfile_input = pysam.AlignmentFile(bam_file_path, "rb")
    output_files_by_barcode = {barcode: pysam.AlignmentFile(output_folder / f"{barcode}.bam", "wb",
                                                            template=samfile_input)
                               for barcode in cell_barcodes}
    for read in tqdm(samfile_input, desc='Processing reads from the input .bam file'):
        if not (read.has_tag('CB') and read.has_tag('UB')):
            continue

        read_barcode = read.get_tag('CB')
        if read_barcode in cell_barcodes:
            output_files_by_barcode[read_barcode].write(read)

    for samfile_output in output_files_by_barcode.values():
        samfile_output.close()
    for barcode in tqdm(cell_barcodes, desc='Indexing .bam files'):
        pysam.index(str(output_folder / f"{barcode}.bam"))


def split_aging_mouse_data_by_cell_barcodes(aging_mouse_data_folder_path: Path) -> None:
    barcodes_df = pd.read_csv(
        aging_mouse_data_folder_path / Path('outs/filtered_gene_bc_matrices/gencode.vM19/barcodes.tsv'),
        delimiter='\t',
        names=['barcode'])
    split_bam_file_by_cell_barcodes(bam_file_path=aging_mouse_data_folder_path / 'outs/possorted_genome_bam.bam',
                                    cell_barcodes=set(barcodes_df['barcode']),
                                    output_folder=split_by_cells_folder_path / aging_mouse_data_folder_path.name)


if __name__ == "__main__":
    split_aging_mouse_data_by_cell_barcodes(raw_data_folder_path / '10X_P5_0')

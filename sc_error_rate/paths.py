from pathlib import Path

project_path = Path(__file__).parent.parent.resolve()
data_folder_path = project_path / 'data'

raw_data_folder_path = data_folder_path / 'raw_data'
split_by_cells_folder_path = data_folder_path / 'split_by_cells'

reference_genome_folder_path = data_folder_path / 'reference_genome'

detected_errors_data_folder_path = data_folder_path / 'detected_errors_data'

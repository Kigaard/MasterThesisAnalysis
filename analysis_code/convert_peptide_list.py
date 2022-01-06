import pandas as pd
import pathlib

import tqdm


def _convert_peptide_list_file(peptide_list_file) -> None:
    """
    Convert peptide list file from the old Excel (XLS) to new Excel (XLSX) format.
    The files are created in the same directory.

    :param peptide_list_file: The peptide list.
    """
    # Create the output name
    output_file_name = peptide_list_file.with_suffix(".xlsx")
    # Load in the file
    peptide_list: pd.DataFrame = pd.read_excel(peptide_list_file, engine="xlrd")
    # Select unique peptides only and set the index to the spectrum id
    peptide_list = peptide_list[peptide_list['V'] == 'Y']
    peptide_list = peptide_list.set_index("spec id")
    # Save the list to an excel file.
    peptide_list.to_excel(output_file_name)


def convert_peptide_list_files(peptide_list_directory: str) -> None:
    """
    Convert peptide list files in a directory from the the old Excel (XLS) to new Excel (XLSX) format.
    The files are created in the same directory.

    :param peptide_list_directory: The peptide list directory.
    """
    # Get the files
    directory = pathlib.Path(peptide_list_directory)
    files = [directory / file for file in directory.glob("*.xls") if file.is_file()]
    # Convert each of the files
    for file in tqdm.tqdm(files):
        _convert_peptide_list_file(file)

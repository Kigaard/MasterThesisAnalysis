import pandas as pd
import os

import tqdm


def convert_peptide_list(peptide_list_file: str):
    """
    Convert peptide list file.

    :param peptide_list_file: The peptide list
    """
    output_file_name = os.path.splitext(peptide_list_file)[0] + ".xlsx"
    peptide_list: pd.DataFrame = pd.read_excel(peptide_list_file, engine="xlrd")
    peptide_list = peptide_list[peptide_list['V'] == 'Y']
    peptide_list = peptide_list.set_index("spec id")
    peptide_list.to_excel(output_file_name)


if __name__ == '__main__':

    directory = r""
    files = [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(".xls")]

    for file in tqdm.tqdm(files):
        convert_peptide_list(file)


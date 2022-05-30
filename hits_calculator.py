import pathlib
from typing import Dict

import pandas as pd


def calculate_hit_statistics(directory: pathlib.Path) -> pd.DataFrame:
    """
    Calculate hit statics

    :param directory: The directory.
    :return: The pandas data frame with the statistics.
    """
    files: list[tuple[str, pathlib.Path]] = [(file.stem, directory / file) for file in directory.glob("*.xlsx")
                                             if file.is_file()]

    result_dict: dict = {}

    for file_name, file_path in files:
        data: pd.DataFrame = pd.read_excel(file_path, usecols=["seq", "#", "modifs"])

        # Check C-terminal modification for trypsin digestions
        if "tryp" in file_name:
            continue
            #o18_1 = data["modifs"].str.contains("@2.004")
            #o18_2 = data["modifs"].str.contains("@4.009")
            #data = data[o18_1 | o18_2]

        result_dict[file_name] = _calculate_file_hit_statistics(hits_df=data)

    return pd.DataFrame.from_dict(result_dict, orient="index")


def _calculate_file_hit_statistics(hits_df: pd.DataFrame) -> dict:
    """
    Calculate the hit statistics for an single file.

    :param hits_df: The file path to the hit file.
    :return: The dictionary with the statistics.
    """
    result_dict: Dict[str, int] = {'TotalHits': 0, 'UniqueHits': 0, 'HitsOnlySeq': 0}

    result_dict["TotalHits"] = sum(hits_df["#"])
    result_dict["UniqueHits"] = len(hits_df.index)
    result_dict["HitsOnlySeq"] = len(list(set(hits_df["seq"])))

    return result_dict


def main():
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)

    file_path: pathlib.Path = \
        pathlib.Path(
            r"C:\Users\spec-makie17\Documents\Experiments\FigureGeneration\Bait_Take3\Lists_CBM")

    df = calculate_hit_statistics(file_path)
    print(df)
    df.to_excel(file_path / "../HitsCBM.xlsx")


if __name__ == "__main__":
    main()

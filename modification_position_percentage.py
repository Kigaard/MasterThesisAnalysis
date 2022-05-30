import functools
import pathlib

import numpy as np
import pandas as pd


def read_filter_dataframes(file_path: pathlib.Path, portion: str) -> tuple[
    pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    def read_filter_df(batch_name: str) -> pd.DataFrame:
        df: pd.DataFrame = pd.read_excel(file_path, sheet_name=f"Batch {batch_name} - {portion}", index_col=0)
        positions = df.index
        df["Position"] = positions
        position_integer = [int(str.replace(i, "M", "").replace("P", "").replace("H", "")) for i in positions]
        df["PositionInt"] = position_integer
        df = df.set_index("PositionInt")

        df = df[df["TotalModSpectra"] > 1]
        del df["TotalModSpectra"]
        del df["ModifiedSpectra"]

        df = df[["Position", "Percentage"]]
        return df

    batch25_df: pd.DataFrame = read_filter_df("2.5")
    batch35_df: pd.DataFrame = read_filter_df("3.5")

    batch1903_df: pd.DataFrame = read_filter_df("1903F")
    batch1806_df: pd.DataFrame = read_filter_df("1806F")

    return batch25_df, batch35_df, batch1903_df, batch1806_df


def combine_dataframes(dfs: tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]):
    new_columns = ["Position", "2.5", "3.5", "19/03F", "18/06F"]
    df_left = functools.reduce(lambda left, right: pd.merge(left, right, how="left", on='Position'), dfs)
    df_left.columns = new_columns
    df_right = functools.reduce(lambda left, right: pd.merge(left, right, how="right", on='Position'), dfs)
    df_right.columns = new_columns

    combined_df = pd.merge(df_left, df_right, how="outer", on="Position")
    combined_df = combined_df.set_index("Position")
    combined_df["2.5"] = np.nan
    combined_df["3.5"] = np.nan
    combined_df["19/03F"] = np.nan
    combined_df["18/06F"] = np.nan

    combined_df = combined_df.fillna(0)

    for _, row in combined_df.iterrows():
        if row["2.5_x"] != 0:
            row["2.5"] = row["2.5_x"]
        elif row["2.5_y"] != 0:
            row["2.5"] = row["2.5_y"]

        if row["3.5_x"] != 0:
            row["3.5"] = row["3.5_x"]
        elif row["3.5_y"] != 0:
            row["3.5"] = row["3.5_y"]

        if row["19/03F_x"] != 0:
            row["19/03F"] = row["19/03F_x"]
        elif row["19/03F_y"] != 0:
            row["19/03F"] = row["19/03F_y"]

        if row["18/06F_x"] != 0:
            row["18/06F"] = row["18/06F_x"]
        elif row["18/06F_y"] != 0:
            row["18/06F"] = row["18/06F_y"]

    combined_df = combined_df.loc[:, combined_df.columns.isin(["2.5", "3.5", "19/03F", "18/06F"])]
    combined_df.replace(0, np.nan, inplace=True)

    combined_df = combined_df[combined_df.iloc[:, 0:3].notna().sum(axis=1) > 1]

    return combined_df


def main():
    portion: str = "Flowthrough"
    file: pathlib.Path = pathlib.Path(
        r"C:\Users\spec-makie17\Documents\Experiments\220425_Batches_PH22006\Oxidations.xlsx")
    dfs = read_filter_dataframes(file_path=file, portion=portion)
    combined = combine_dataframes(dfs)
    combined["Average"] = round(combined.mean(axis=1), 3)

    combined.to_excel(file.parent / f"Cys_Combined_{portion}.xlsx")


if __name__ == "__main__":
    main()

import ast
import re
import pandas as pd
import numpy as np
from pathlib import Path


def parse_lattice(filepath) -> pd.DataFrame:
    return pd.read_csv(filepath, header=None)


def parse_food_data(datapath, time_filter: list = None) -> pd.DataFrame:
    df = _parse_dfs(datapath, time_filter)
    df["sigma_list"] = df["sigma_list"].apply(lambda string: [int(x) for x in string.split(' ')])
    return df


def parse_cell_data(datapath,
                    time_filter: list = None) -> pd.DataFrame:
    """Parses a CSV file into a dataframe containing cell information.

    The format of this dataframe is essential to API stability, please be mindful of changes to it.
    """
    celldf = _parse_dfs(datapath, time_filter)
    list_attrs = ["neighbour_list", "neighbourJ_list", "in_scale_list", "reg_threshold_list",
                  "reg_w_innode_list", "reg_w_regnode_list", "out_threshold_list",
                  "out_w_regnode_list"]
    for key in list_attrs:
        # If the string only contains one object it will be already parsed as a literal type
        # so no need to split it
        if celldf[key].dtype == object:
            celldf[key] = celldf[key].apply(
                lambda string: [ast.literal_eval(x) for x in string.split(' ')]
            )
    celldf = celldf.sort_values(["time", "sigma"])

    if Path(datapath).is_file():
        return celldf.set_index("sigma", drop=False)
    # Multi-index
    return celldf.set_index(["time", "sigma"], drop=False)


def get_time_points(datapath):
    times = set()
    for filepath in Path(datapath).iterdir():
        m = re.search(r"t(\d+).csv", filepath.name)
        if m is not None:
            times.add(int(m.group(1)))
    return sorted(times)


def build_time_filter(time_points, start=0, end=float("inf"), n=None):
    """Selects 'n' interspaced time points in the range ['start', 'end']."""
    time_points = np.unique(time_points)
    time_points = time_points[(time_points >= start) & (time_points <= end)]

    if n is None:
        n = len(time_points)
    else:
        n = min(n, len(time_points))

    indexes = np.round(np.linspace(0, len(time_points) - 1, n)).astype(int)
    return time_points[indexes]


def _parse_dfs(datapath, time_filter):
    datapath = Path(datapath)

    if datapath.is_file():
        return pd.read_csv(datapath).sort_values("time")

    dfs = []
    for filepath in datapath.iterdir():
        df = pd.read_csv(filepath)
        if time_filter is None or df["time"][0] in time_filter:
            dfs.append(df)
    return pd.concat(dfs).sort_values("time")
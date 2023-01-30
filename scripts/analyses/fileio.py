import re
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from enlighten import Counter
from make_netgraphs import _gene_attrs

logger = logging.getLogger(__name__)
_str_attrs = _gene_attrs + ["neighbour_list", "Jneighbour_list"]


def parse_lattice(filepath) -> pd.DataFrame:
    logger.info(f"Parsing dataframe from: '{filepath}'")
    return pd.read_csv(filepath, header=None)


def parse_food_data(datapath, time_filter: list = None, trust_filenames=True) -> pd.DataFrame:
    df = _parse_dfs(datapath, time_filter, trust_filenames)
    return df.sort_values("time")


def parse_cell_data(datapath,
                    time_filter: list = None,
                    trust_filenames=True) -> pd.DataFrame:
    """Parses a CSV file into a dataframe containing cell information.

    The format of this dataframe is essential to API stability, please be mindful of changes to it.
    
    :param time_filter:
    :param datapath:
    :param trust_filenames: Whether to speed-up the filtering process by trusting that the file
    names represent the time of the simulation.
    """

    celldf = _parse_dfs(datapath,
                        time_filter,
                        trust_filenames,
                        dtype={col: str for col in _str_attrs})
    celldf = celldf.sort_values(["time", "sigma"])

    # Can't drop the columns because they might be needed by other functions!
    if Path(datapath).is_file():
        celldf = celldf.set_index("sigma", drop=False)
        # Some operations dont like indexes and columns with the same name
        celldf.index.name = "sigma_i"
        return celldf
    # Multi-index
    celldf = celldf.set_index(["time", "sigma"], drop=False)
    celldf.index.names = ["time_i", "sigma_i"]
    return celldf


def parse_grave_data(datapath, time_filter=None, trust_filenames=True) -> pd.DataFrame:
    gravedf = _parse_dfs(datapath, time_filter, trust_filenames)
    gravedf = gravedf.sort_values(["time_death", "sigma"])
    return gravedf.set_index(["time_death", "sigma"], drop=False)


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


def write_plot(fig, outputfile):
    logger.info(f"Writing plot to: {outputfile}")
    if outputfile[-5:] == ".html":
        fig.write_html(outputfile)
    else:
        fig.write_image(outputfile)


def _parse_dfs(datapath, time_filter, trust_filenames, **csv_kwargs):
    logger.info(f"Parsing dataframe(s) from: '{datapath}'")    
    datapath = Path(datapath)

    if datapath.is_file():
        return pd.read_csv(datapath)

    dfs = []
    filepaths = list(datapath.iterdir())
    pbar = Counter(total=len(filepaths), desc="CSV files iterated")
    for filepath in filepaths:
        if time_filter is None:
            dfs.append(pd.read_csv(filepath, **csv_kwargs))
        # Skips unecessary calls to pd.read_csv by trusting file name
        elif not trust_filenames:
            df = pd.read_csv(filepath, **csv_kwargs)
            if df["time"][0] in time_filter:
                dfs.append(df)
        elif int(re.search(r"\d+", filepath.name).group()) in time_filter:
            dfs.append(pd.read_csv(filepath, **csv_kwargs))
        pbar.update()
    return pd.concat(dfs)
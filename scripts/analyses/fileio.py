import re
import logging
import pandas as pd
import numpy as np
import plotly.express as px
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
    total = len(time_filter) if time_filter is not None else len(filepaths)
    pbar = Counter(total=total, desc="CSV files iterated")
    for filepath in filepaths:
        # Some very confusing conditionals to speed up processing
        # Dont mess with this, it was very hard to figure out the optimal way to do it
        add = False
        if time_filter is None:
            add = True
        elif trust_filenames:
            if int(re.search(r"\d+", filepath.name).group()) not in time_filter:
                continue
            add = True
        df = pd.read_csv(filepath, **csv_kwargs)
        if add or df["time"][0] in time_filter:
            dfs.append(df)
            pbar.update()
    return pd.concat(dfs)


def plot_lattice(latdf, outputfile, colors=None, bg_color=None):
    """Plots the lattice as an image.
    :param outputfile:
    :param bg_color:
    :param colors: If given, will be used to map each sigma to a value on the image. Should be a
        dictionary or pd.Series with sigmas as indexes and integers or rgb values as values.
    :param latdf:
    """
    lat = latdf.values
    if colors is None:
        sigs = np.unique(lat)
        # Remove zeros
        colors = pd.Series(np.random.random(sigs.size - 1), sigs[1:])
    else:
        colors = pd.Series(colors)

    if bg_color is None:
        if len(colors.shape) == 2:
            bg_color = np.zeros(colors.shape[1])
        else:
            bg_color = 0

    colors = colors.reindex(np.arange(colors.index.max() + 1), fill_value=bg_color)
    cvals = colors.values

    img = cvals[lat]
    fig = px.imshow(img)
    write_plot(fig, outputfile)
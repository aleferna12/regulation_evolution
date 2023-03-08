import logging
import argparse
import subprocess
import ast
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
from tempfile import NamedTemporaryFile
from plotly.subplots import make_subplots
from enlighten import Counter
from colorir import *
from scripts.fileio import write_plot, parse_cell_data
from scripts import PROJECT_DIR

logger = logging.getLogger(__name__)


def get_parser():
    def run(args):
        celldfs = []
        column_titles = []
        for filepath in Path(args.datadir).iterdir():
            celldfs.append(parse_cell_data(filepath))
            column_titles.append(filepath.name.replace(".csv", ""))
        genomesdf, keylock_tups = sample_genomes_and_keylocks(celldfs, args.cell_index)
        sweepdf = sweep_genomes(genomesdf,
                                args.min_chem,
                                args.max_chem,
                                args.step_chem,
                                args.min_foodparc,
                                args.max_foodparc,
                                args.step_foodparc,
                                args.mcss)
        
        marker_style = {}
        if args.marker_style:
            for style_name, style_val in zip(args.marker_style[::2], args.marker_style[1::2]):
                try:
                    marker_style[style_name] = ast.literal_eval(style_val)
                except ValueError:  # Assume its a string
                    marker_style[style_name] = style_val
        fig = plot_genome(sweepdf,
                          args.outputfile,
                          args.Jmed,
                          args.Jalpha,
                          args.weights,
                          keylock_tups,
                          column_titles=column_titles,
                          two_colorscales=args.two_colorscales,
                          colorscale=args.colorscale,
                          sep=args.sep,
                          marker_spacing=args.marker_spacing,
                          mig_marker_color=args.mig_marker_color,
                          div_marker_color=args.div_marker_color,
                          **marker_style)
        if args.show:
            fig.show("browser")

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Performs parameter sweeps in genomes sampled from simulations"
    )
    parser.add_argument(
        "datadir",
        help="Directory containing the cell CSV files. "
             "A random genome will be sampled from each of the files in "
             "this directory and queried with a range of inputs. You can select "
             "which genome will be queried by having a single row per file"
    )
    parser.add_argument(
        "outputfile",
        help="Output HTML or SVG file. In the first row, "
             "the gamma of each parameter coordinate is plotted by comparing "
             "its own key and lock. In the second row, the key and lock of each parameter "
             "coordinate is compared to the most common keys and locks among cells with "
             "the alternate cell type in each of the simulations. "
             "Column titles are assigned according to the names of the "
             "cell data files in 'datadir'"
    )
    parser.add_argument(
        "--Jmed",
        default=14,
        type=int,
        help="J value against the medium used in the gamma calculations"
    )
    parser.add_argument(
        "--Jalpha",
        default=7,
        type=int,
        help="J value which is summed to cell adhesion strengths in the gamma calculations"
    )
    parser.add_argument(
        "-w",
        "--weights",
        nargs="+",
        default=list(range(1, 7)),
        type=float,
        help="Weights applied to the keys and locks for gamma calculation"
    )
    parser.add_argument(
        "--min-chem",
        default=40,
        type=float,
        help="Minimum chemotactical signal used as input in the genome sweep"
    )
    parser.add_argument(
        "--max-chem",
        default=90,
        type=float,
        help="Maximum chemotactical signal used as input in the genome sweep"
    )
    parser.add_argument(
        "--step-chem",
        default=1,
        type=float,
        help="Step with which to increase the chem signal "
             "input for the genome weep"
    )
    parser.add_argument(
        "--min-foodparc",
        default=0,
        type=float,
        help="Minimum number of food parcels used as input in the genome sweep"
    )
    parser.add_argument(
        "--max-foodparc",
        default=25,
        type=float,
        help="Maximum number of food parcels used as input in the genome sweep"
    )
    parser.add_argument(
        "--step-foodparc",
        default=0.5,
        type=float,
        help="Step with which to increase the number of food parcels "
             "input for the genome weep"
    )
    parser.add_argument(
        "--mcss",
        default=50,
        type=int,
        help="How many times update the genome with the same inputs. "
             "Setting this parameter to a very low value causes artifacts in the plots "
             "(less than 10 is kinda dangerous, less than 5 is unreliable)."
    )
    parser.add_argument(
        "--cell-index",
        help="Can be used to chose which cell's genome to plot from each cel data frame (instead of picking one at random). "
             "All data frames must have this index",
        type=int
    )
    parser.add_argument("-s",
                        "--show",
                        action="store_true",
                        help="Show figure after plotting")
    parser.add_argument(
        "--two-colorscales",
        action="store_true",
        help="Use different colorscales for migrating and dividing cells. Otherwise, 'colorscale' will be used for both"
    )
    parser.add_argument(
        "--colorscale",
        default="Viridis",
        help="If not 'two-colorscales', then which colorscale to use. Any plotly colorscales are available"
    )
    parser.add_argument(
        "--sep",
        default=0,
        type=float,
        help="Gap between cells in the heatmaps"
    )
    parser.add_argument(
        "--marker-spacing",
        default=5,
        type=int,
        help="How much to space the markers apart (in terms of rows/columns in the graph, not space!)."
             "Set to 0 to disable markers"
    )
    parser.add_argument(
        "--mig-marker-color",
        default="#dd4f33",
        help="Color used for the migrating markers"
    )
    parser.add_argument(
        "--div-marker-color",
        default="#4675f0",
        help="Color used for the migrating markers"
    )
    parser.add_argument(
        "--marker-style",
        nargs="+",
        help="Styling properties for the markers passed as pairs <property> <value>. "
             "Look at the docs of plotly's Scatter for a full API. "
             "E.g.: --marker-style marker_symbol square marker_size 5"
    )
    parser.set_defaults(run=run)
    return parser


def sweep_genomes(genomesdf: pd.DataFrame,
                  min_chem,
                  max_chem,
                  step_chem,
                  min_foodparc,
                  max_foodparc,
                  step_foodparc,
                  mccs: int,
                  sweep_bin: str = None):
    logger.info("Starting genome sweep subprocess")
    if sweep_bin is None:
        sweep_bin = PROJECT_DIR / "bin" / "sweep_genomes"
    else:
        sweep_bin = Path(sweep_bin).resolve()
    if not sweep_bin.is_file():
        raise FileNotFoundError("could not resolve path to sweep_genomes binary file, "
                                "build this target with cmake (or specify it's location if "
                                "it's already built)")

    with NamedTemporaryFile() as inputfile, NamedTemporaryFile("r") as outputfile:
        genomesdf.to_csv(inputfile, index=False)
        proc = subprocess.run([str(arg) for arg in [sweep_bin,
                                                    inputfile.name,
                                                    outputfile.name,
                                                    min_chem,
                                                    max_chem,
                                                    step_chem,
                                                    min_foodparc,
                                                    max_foodparc,
                                                    step_foodparc,
                                                    mccs]], capture_output=True)
        if proc.returncode:
            raise RuntimeError(proc.stderr)
        return pd.read_csv(outputfile)


def sample_genomes_and_keylocks(celldfs: "list[pd.DataFrame]", cell_index=None):
    logger.info("Sampling information from cell dataframes")
    genomes = []
    keylock_tups = []
    for celldf in celldfs:
        if cell_index is None:
            genomes.append(celldf.sample(1))
        else:
            genomes.append(celldf.iloc[[cell_index]])
        # If two values are found in the same freq, pick the first one
        kldf = celldf.groupby("tau")[["jkey_dec", "jlock_dec"]].agg(lambda s: pd.Series.mode(s)[0])
        keylock_tups.append(tuple(kldf.values.flatten()))
    genomes = pd.concat(genomes)[["innr",
                                  "regnr",
                                  "outnr",
                                  "in_scale_list",
                                  "reg_threshold_list",
                                  "reg_w_innode_list",
                                  "reg_w_regnode_list",
                                  "out_threshold_list",
                                  "out_w_regnode_list"]]
    return genomes, keylock_tups


def plot_genome(sweepdf: pd.DataFrame,
                outputfile,
                Jmed,
                Jalpha,
                weights,
                keylock_tups: "list[tuple[int, int, int, int]]",
                column_titles,
                two_colorscales,
                colorscale,
                sep,
                marker_spacing,
                mig_marker_color,
                div_marker_color,
                **marker_style):
    logger.info(f"Plotting the network outputs to: {outputfile}")
    weights = np.array(weights, dtype=float)
    gdf = sweepdf.groupby("id")
    if len(gdf) != len(keylock_tups):
        raise ValueError("received a different number of genomes and keys/locks values "
                         f"(respectively {len(gdf)} and {len(keylock_tups)})")

    fig = make_subplots(
        2,
        len(gdf),
        shared_xaxes=True,
        column_titles=column_titles
    )
    pbar = Counter(total=len(gdf), desc="Genomes analyzed")
    for (g_id, df), keylocks in zip(gdf, keylock_tups):
        traces = make_heatmaps(df,
                               Jmed,
                               Jalpha,
                               weights,
                               keylocks,
                               two_colorscales,
                               colorscale,
                               sep,
                               marker_spacing,
                               mig_marker_color,
                               div_marker_color,
                               **marker_style)
        fig.add_traces(traces[:4], rows=[1, 1, 2, 2], cols=g_id + 1)
        fig.add_traces(traces[4:], rows=1, cols=g_id + 1)
        fig.add_traces(traces[4:], rows=2, cols=g_id + 1)
        fig.update_xaxes(range=[traces[0].x[0], traces[0].x[-1]],
                         constrain="domain",
                         showticklabels=True)
        fig.update_yaxes(scaleanchor="x" + ("" if not g_id else str(g_id + 1)),
                         scaleratio=.5,
                         constrain="domain",
                         range=[traces[0].y[0], traces[0].y[-1]],
                         col=g_id + 1)
        pbar.update()
    fig.update_xaxes(title="food parcels", row=2)
    fig.update_yaxes(title="chemotactic signal", col=1)
    fig.update_layout(width=400 * (len(gdf)) + two_colorscales * 100, height=600)
    write_plot(fig, outputfile)
    return fig


def make_heatmaps(df, 
                  Jmed,
                  Jalpha, 
                  weights,
                  keylocks: "tuple[int, int, int, int]",
                  two_colorscales,
                  colorscale,
                  sep,
                  marker_spacing,
                  mig_marker_color,
                  div_marker_color,
                  **marker_style):
    chem_vals = df.chem.unique()
    foodparc_vals = df.food.unique()
    img_size = (len(chem_vals), len(foodparc_vals))

    max_gamma = get_gamma(
        np.ones(weights.size),
        np.zeros(weights.size),
        np.ones(weights.size),
        np.zeros(weights.size),
        Jmed,
        Jalpha,
        weights
    )
    min_gamma = get_gamma(
        np.zeros(weights.size),
        np.zeros(weights.size),
        np.zeros(weights.size),
        np.zeros(weights.size),
        Jmed,
        Jalpha,
        weights
    )
    jkeys = [bistring_from_dec(jkey_dec, weights.size)
             for jkey_dec in df.jkey_dec]
    jlocks = [bistring_from_dec(jlock_dec, weights.size)
              for jlock_dec in df.jlock_dec]
    self_gammas = [get_gamma(k, l, k, l, Jmed, Jalpha, weights) for k, l in zip(jkeys, jlocks)]
    sel_gammas_img = np.reshape(self_gammas, img_size)
    tau_img = df.tau.values.reshape(img_size)

    # Handpicked colors, dont mess with it for gods sake
    # If it becomes hard to see if unicellular cells are mig or div, replace the first color
    # of each grad with #ffeee6 and #edf5fc repectively
    migscale = PolarGrad(["#f9f0e3", "#e8ba4c", "#a43829", "#59000b"],
                         color_coords=[0, .5, .75, 1]).to_plotly_colorscale()
    divscale = PolarGrad(["#f6f2fa", "#bf94ec", "#386388", "#001f4d"],
                         color_coords=[0, .5, .75, 1]).to_plotly_colorscale()
    mig_img = np.where(tau_img == 1, sel_gammas_img, np.nan)
    div_img = np.where(tau_img == 2, sel_gammas_img, np.nan)

    traces = []
    # Self-gamma
    traces.append(go.Heatmap(
        z=mig_img,
        x=foodparc_vals,
        y=chem_vals,
        zmin=min_gamma,
        zmax=max_gamma,
        xgap=sep,
        ygap=sep,
        colorbar=go.heatmap.ColorBar(y=1, yanchor="top", len=0.5),
        colorscale=migscale if two_colorscales else colorscale
    ))
    traces.append(go.Heatmap(
        z=div_img,
        x=foodparc_vals,
        y=chem_vals,
        zmin=min_gamma,
        zmax=max_gamma,
        xgap=sep,
        ygap=sep,
        colorbar=go.heatmap.ColorBar(xpad=75, y=1, yanchor="top", len=0.5),
        colorscale=divscale if two_colorscales else colorscale,
        showscale=two_colorscales,
    ))

    migkey = bistring_from_dec(keylocks[0], weights.size)
    miglock = bistring_from_dec(keylocks[1], weights.size)
    divkey = bistring_from_dec(keylocks[2], weights.size)
    divlock = bistring_from_dec(keylocks[3], weights.size)
    other_gammas = []
    for tau, jkey, jlock in zip(df.tau, jkeys, jlocks):
        if tau == 1:
            other_key = divkey
            other_lock = divlock
        else:
            other_key = migkey
            other_lock = miglock
        other_gammas.append(get_gamma(jkey, jlock, other_key, other_lock, Jmed, Jalpha, weights))
    other_gammas_img = np.reshape(other_gammas, img_size)
    traces.append(go.Heatmap(
        z=np.where(tau_img == 1, other_gammas_img, np.nan),
        x=foodparc_vals,
        y=chem_vals,
        zmin=min_gamma,
        zmax=max_gamma,
        xgap=sep,
        ygap=sep,
        colorscale=migscale if two_colorscales else colorscale,
        showscale=False
    ))
    traces.append(go.Heatmap(
        z=np.where(tau_img == 2, other_gammas_img, np.nan),
        x=foodparc_vals,
        y=chem_vals,
        zmin=min_gamma,
        zmax=max_gamma,
        xgap=sep,
        ygap=sep,
        colorscale=divscale if two_colorscales else colorscale,
        showscale=False
    ))

    if marker_spacing:
        symbol_img = np.zeros_like(tau_img)
        symbol_img[::marker_spacing, ::marker_spacing] = tau_img[::marker_spacing, ::marker_spacing]
        symbol_img = symbol_img.T
        mig_indexes = np.nonzero(symbol_img == 1)
        div_indexes = np.nonzero(symbol_img == 2)
        traces.append(go.Scatter(
            x=foodparc_vals[mig_indexes[0]],
            y=chem_vals[mig_indexes[1]],
            mode="markers",
            showlegend=False,
            marker_color=mig_marker_color,
            marker_line_color=mig_marker_color,
            **marker_style
        ))
        traces.append(go.Scatter(
            x=foodparc_vals[div_indexes[0]],
            y=chem_vals[div_indexes[1]],
            mode="markers",
            showlegend=False,
            marker_color=div_marker_color,
            marker_line_color=div_marker_color,
            **marker_style
        ))
    return traces


def get_median_gamma(key_locks, Jmed, Ja, weights, n=None):
    """Calculates the medium gamma from all x all cells.

    :param key_locks: Numpy matrix of (key, lock) pairs or each cell
    :param n: Samples n key-lock pairs from total population to speed up calculations
    :param Ja:
    :param Jmed:
    :param weights:
    """
    if n is not None:
        rng = np.random.default_rng(0)
        indexes = rng.choice(len(key_locks), min(n, len(key_locks)), replace=False)
        key_locks = key_locks[indexes]

    gammas = []
    for k1, l1 in key_locks:
        for k2, l2 in key_locks:
            gammas.append(get_gamma(k1, l1, k2, l2, Jmed, Ja, weights))

    return np.median(gammas)


def get_cell_cell_j(k1, l1, k2, l2, Ja, weights):
    return Ja + np.sum((k1 == l2) * weights) + np.sum((k2 == l1) * weights)


def get_gamma(k1, l1, k2, l2, Jmed, Ja, weights):
    return Jmed - get_cell_cell_j(k1, l1, k2, l2, Ja, weights) / 2


def bistring_from_dec(dec, width):
    return np.fromiter(np.binary_repr(dec, width=width), dtype=int)


def test():
    parser = get_parser()
    args = parser.parse_args(["scripts/genome/data", "scripts/genome/plots.svg", "-s"])
    args.run(args)


if __name__ == "__main__":
    test()

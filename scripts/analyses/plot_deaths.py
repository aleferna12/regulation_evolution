import argparse
import logging
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from fileio import *

logger = logging.getLogger(__name__)


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(prog="plot_deaths",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("datadir", help="directory containing the cell grave data frames")
    parser.add_argument("outputfile", help="SVG output file")
    parser.add_argument("-b",
                        "--bin-size",
                        help="determine the plot bin size in MCS",
                        type=int,
                        default=50000)
    args = parser.parse_args()

    gravedf = parse_grave_data(args.datadir)
    plot_deaths(gravedf, args.outputfile, args.bin_size)

    logger.info("Finished")


def plot_deaths(gravedf: pd.DataFrame, outputfile, bin_size):
    bins = range(0, gravedf["time_death"].max() + bin_size, bin_size)
    intervals = pd.cut(gravedf["time_death"], bins=bins)
    deaths_interval = gravedf.groupby(intervals).size()
    x = ["%i - %i kMCS" % (ti.left / 1000, ti.right / 1000)
         for ti in intervals.cat.categories]

    fig = make_subplots(3,
                        1,
                        shared_xaxes=True,
                        subplot_titles=["Number of deaths",
                                        "Cause of death in %",
                                        "Mean age of death"])

    # Plot reason of death
    reasons = gravedf.groupby(["reason", intervals]).size()
    for i, reason in enumerate(gravedf["reason"].unique()):
        color = px.colors.qualitative.Plotly[i]
        fig.add_trace(go.Scatter(
            x=x,
            y=reasons[reason],
            stackgroup=0,
            line_color=color,
            name=reason
        ), row=1, col=1)
        fig.add_trace(go.Scatter(
            x=x,
            y=reasons[reason] / deaths_interval,
            stackgroup=1,
            line_color=color,
            showlegend=False
        ), row=2, col=1)

    # Plot average age of death
    fig.add_trace(go.Scatter(
        x=x,
        y=gravedf["time_since_birth"].groupby(intervals).mean(),
        name="age",
        showlegend=False,
    ), row=3, col=1)

    fig.update_yaxes(range=[0, 1], row=2)
    fig.update_xaxes(range=[0, len(x) - 1])

    write_plot(fig, outputfile)


if __name__ == "__main__":
    main()

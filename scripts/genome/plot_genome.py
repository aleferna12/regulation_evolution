import logging
import argparse

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from plotly.subplots import make_subplots
from scripts.fileio import write_plot
from scripts.analyses.plot_timeline import get_gamma


def get_parser():
    def run():
        pass

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "inputfile",
        help="Input file with the outputs of a genome. Generate these with 'sweep_genome.cpp'."
    )
    parser.add_argument("outputfile", help="Output HTML or SVG file.")
    parser.set_defaults(run=run)
    return parser


def plot_genome(inputfile, outputfile):
    SEP = 0
    JMED = 14
    JALPHA = 7
    WEIGHTS = np.arange(1, 7)

    df = pd.read_csv(inputfile)
    _, adh_strats = np.unique(np.stack([df.jkey_dec, df.jlock_dec], axis=1),
                              axis=0,
                              return_inverse=True)

    chem_vals = df.chem.unique()
    foodparc_vals = df.food.unique()
    img_size = (len(chem_vals), len(foodparc_vals))

    jkeys = [np.fromiter(np.binary_repr(jkey_dec, width=WEIGHTS.size), dtype=int)
             for jkey_dec in df.jkey_dec]
    jlocks = [np.fromiter(np.binary_repr(jlock_dec, width=WEIGHTS.size), dtype=int)
              for jlock_dec in df.jlock_dec]
    gammas = [get_gamma(k, l, k, l, JMED, JALPHA, WEIGHTS) for k, l in zip(jkeys, jlocks)]
    max_gamma = get_gamma(
        np.ones(WEIGHTS.size),
        np.zeros(WEIGHTS.size),
        np.ones(WEIGHTS.size),
        np.zeros(WEIGHTS.size),
        JMED,
        JALPHA,
        WEIGHTS
    )
    min_gamma = get_gamma(
        np.zeros(WEIGHTS.size),
        np.zeros(WEIGHTS.size),
        np.zeros(WEIGHTS.size),
        np.zeros(WEIGHTS.size),
        JMED,
        JALPHA,
        WEIGHTS
    )
    gammas_img = np.reshape(gammas, img_size)
    tau_img = df.tau.values.reshape(img_size)
    gammas_tau1 = np.where(tau_img == 1, gammas_img, np.nan)
    gammas_tau2 = np.where(tau_img == 2, gammas_img, np.nan)

    fig = go.Figure()
    fig.add_trace(go.Heatmap(
        z=gammas_tau1,
        x=foodparc_vals,
        y=chem_vals,
        zmin=min_gamma,
        zmax=max_gamma,
        colorscale="Reds",
        colorbar=go.heatmap.ColorBar(x=1, xanchor="right"),
        xgap=SEP,
        ygap=SEP,
    ))
    fig.add_trace(go.Heatmap(
        z=gammas_tau2,
        x=foodparc_vals,
        y=chem_vals,
        zmin=min_gamma,
        zmax=max_gamma,
        colorscale="Blues",
        xgap=SEP,
        ygap=SEP,
    ))
    fig.update_xaxes(constrain="domain",
                     title="food parcels")
    fig.update_yaxes(scaleanchor="x",
                     scaleratio=1,
                     title="chemotactic signal")
    fig.show("browser")


def test():
    plot_genome("data/uniprop_out.csv", "")


if __name__ == "__main__":
    test()
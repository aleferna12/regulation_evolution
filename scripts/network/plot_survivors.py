"""This is an alternative way to plot a genealogy tree."""

import plotly.graph_objects as go
from ete3 import Tree
from colorir import *

_cs = Palette.load()


def main():
    tree = Tree("./data/genealogy_test/trees/tree2.newick")
    plot_survivors(tree)


def plot_survivors(tree, marker_size=5):
    # Add time to root here since it can't be saved in NHX
    tree.add_feature("time", 0)

    season_size = {}
    for node in tree.traverse():
        if node.time not in season_size:
            season_size[node.time] = 0
        season_size[node.time] += 1
        node.add_feature("index", season_size[node.time])
    survivors = max(season_size.values())

    xs = []
    ys = []
    # Reiterate tree to set plot positions
    for node in tree.traverse(strategy="levelorder"):
        x = int(node.time)
        step = (survivors + 1) / (season_size[node.time] + 1)
        y = step * node.index
        # Next gen will need this
        node.add_feature("pos", (x, y))

        if not node.is_root():
            anc_x, anc_y = node.up.pos
            xs += [x, anc_x, None]
            ys += [y, anc_y, None]
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=xs,
        y=ys,
        mode="markers+lines",
        marker_size=marker_size,
        line_width=1,
        line_color=_cs.black,
        marker_color=_cs.blue
    ))
    fig.update_layout(height=2 * marker_size * survivors)
    fig.show()


if __name__ == "__main__":
    main()
"""This is an alternative way to plot a genealogy tree."""

import plotly.express as px
import plotly.graph_objects as go
from typing import List
from plotly.subplots import make_subplots
from ete3 import Tree


def main():
    tree = Tree("./data/genealogy_test/trees/tree2.newick")
    plot_survivors([tree])


def plot_survivors(trees: List[Tree]):
    fig = make_subplots(len(trees), 1)
    for i, tree in enumerate(trees, start=1):
        # Add time to root here since it can't be saved in NHX
        tree.add_feature("time", 0)

        season_size = {}
        for node in tree.traverse():
            if node.time not in season_size:
                season_size[node.time] = 0
            season_size[node.time] += 1

        # Reiterate tree to set plot positions
        for node in tree.traverse(strategy="levelorder"):
            node.add_feature("pos", (int(node.time), season_size[node.time]))
            season_size[node.time] -= 1

        xs = []
        ys = []
        for node in tree.iter_descendants():
            x, y = node.pos
            anc_x, anc_y = node.up.pos
            xs += [x, anc_x, None]
            ys += [y, anc_y, None]

        trace = go.Scatter(
            x=xs,
            y=ys,
            mode="markers+lines",
            marker_size=10,
            line_width=1,
            line_color="#000000",
            marker_color="#0000ff"
        )
        fig.add_trace(trace, row=i, col=1)
    fig.update_layout(height=2500)
    fig.show()


if __name__ == "__main__":
    main()
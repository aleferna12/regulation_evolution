import re
from pathlib import Path
from ete3 import Tree


def read_ancestry(path, dead_ends=True):
    path = Path(path)
    seasons = {}
    for filepath in path.glob("anc_*.txt"):
        cells = []
        with open(filepath) as file:
            lines = file.read().split("\n")
        for line in lines:
            if not line:
                continue
            cell_sigma, anc_sigma = line.split(" ")
            if cell_sigma == "0":
                continue
            cells.append((int(cell_sigma), int(anc_sigma)))
        season_i = int(re.search(r"(?<=t)\d+", filepath.name).group())
        seasons[season_i] = cells
    return make_trees(seasons, dead_ends)


# Remember that although the sigmas are added as names, the same id can refer to DIFFERENT cells
def make_trees(seasons, dead_ends=True, names=True):
    # Order by time backwards
    seasons = {k: seasons[k] for k in sorted(seasons, reverse=True)}
    # This dict holds the descendents of the current season of cells
    prev_anc_child = {}
    for season, cell_batch in seasons.items():
        # This dict holds the ancestors of the current season of cells
        next_anc_child = {}
        for cell_sigma, anc_sigma in cell_batch:
            # Keep dead ends from being added to the trees
            if not dead_ends and prev_anc_child and cell_sigma not in prev_anc_child:
                continue
            node = Tree(name=str(cell_sigma) if names else "")
            node.add_feature("season", season)
            if cell_sigma in prev_anc_child:
                node.children = prev_anc_child[cell_sigma]
            if anc_sigma not in next_anc_child:
                next_anc_child[anc_sigma] = []
            next_anc_child[anc_sigma].append(node)
        prev_anc_child = dict(next_anc_child)

    # Construct root trees
    trees = []
    for anc_sigma, children in prev_anc_child.items():
        root = Tree(name=anc_sigma if names else "")
        # The roots date from the start of the first season
        # This is represented as 0 as not to get confused with the end of the first season
        root.add_feature("season", 0)
        root.children = children
        trees.append(root)
    return trees

import argparse
import logging
import re
from scripts.analyses import make_netgraphs, make_trees, plot_deaths, plot_netgraphs, \
    plot_timeline, plot_tree
from scripts.competition import make_competition, make_templates
from scripts.colortable import make_colortable, print_colortable


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        prog="run_script",
        description="To see help about each script run 'python run_script <script> --help'"
    )
    subparsers = parser.add_subparsers(
        title="scripts",
        required=True,
        description="List of scripts that can be executed via this CLI app"
    )
    for mod in [make_netgraphs,
                make_trees,
                plot_deaths,
                plot_netgraphs,
                plot_timeline,
                plot_tree,
                make_competition,
                make_templates,
                make_colortable,
                print_colortable]:
        subparser = mod.get_parser()
        subparsers.add_parser(
            name=re.search(r"\w+$", mod.__name__).group(),
            parents=[subparser],
            help=subparser.description,
            add_help=False
        )

    args = parser.parse_args()
    args.run(args)


if __name__ == "__main__":
    main()
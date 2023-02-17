import argparse
from colorir import *


def get_parser():
    def run(args):
        pal = read_colortable(args.filepath)
        swatch(pal)

    parser = argparse.ArgumentParser(description="Print a color table file to the terminal")
    parser.add_argument("filepath")
    parser.set_defaults(run=run)
    return parser


def read_colortable(filepath) -> StackPalette:
    with open(filepath) as file:
        pal_raw = file.read()
    pal = StackPalette()
    color_fmt = ColorFormat(sRGB)
    for row in pal_raw.split("\n"):
        if not row:
            continue
        color = tuple(int(x) for x in row.split()[1:])
        pal.add(color_fmt.format(color))
    return pal
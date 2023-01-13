from pathlib import Path
from colorir import *

config.DEFAULT_PALETTES_DIR = Path(__file__).resolve().parent / "palettes"


def main():
    pal = StackPalette.load("categ8") \
          + StackPalette.load("chemgrad64_3") \
          + StackPalette.load("divgrad8") \
          + StackPalette.load("miggrad8")
    write_colortable(
        pal,
        Path(__file__).parent.parent.parent / "data" / "colortable.ctb"
    )


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


def write_colortable(spalette: StackPalette, outfile):
    spalette.color_format = ColorFormat(sRGB, round_to=0)
    out_str = ""
    for i, color in enumerate(spalette):
        out_str += f"{i} {color}\n"
    with open(outfile, 'w') as file:
        file.write(out_str.replace('(', '').replace(')', '').replace(',', ''))


def resize_color_list(colors, ncolors):
    return Grad(colors, color_sys=CIELuv).n_colors(ncolors)


if __name__ == "__main__":
    main()
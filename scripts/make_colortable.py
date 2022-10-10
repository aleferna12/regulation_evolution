from colorir import *

START_COLORS = "000000", "ffffff"
FOOD_COLORS = Palette.load("rainbow").colors[::-1] + ["000000"]
END_COLORS = "ffffff", "000000"
OUTFILE = "../data/sunset.ctb"


def make_ncolors(colors, ncolors):
    grad = Grad(colors,
                color_sys=CIELuv,
                color_format=ColorFormat(sRGB, round_to=0))
    return grad.n_colors(ncolors)


def main():
    start_vals = [""] + make_ncolors(START_COLORS, 14)
    food_vals = make_ncolors(FOOD_COLORS, 29)
    end_vals = make_ncolors(END_COLORS, 211)
    file_str = ""
    for i, color in enumerate(start_vals + food_vals + end_vals):
        file_str += f"{i} {' '.join(str(val) for val in color)}\n"
    with open(OUTFILE, 'w') as file:
        file.write(file_str)


if __name__ == "__main__":
    main()
from colorir import *
import plotly.express as px

NCOLORS = 50
COLORS = "#fefe00", "#fe0000"
OUTFILE = "../data/yl_rd.ctb"

grad = Grad(COLORS,
            color_sys=CIELuv,
            color_format=ColorFormat(sRGB, round_to=0))

file_str = ""
for i, color in enumerate(grad.n_colors(NCOLORS)):
    file_str += f"{i} {' '.join(str(val) for val in color)}\n"
with open(OUTFILE, 'w') as file:
    file.write(file_str)
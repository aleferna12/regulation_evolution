import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from colorir import *

config.DEFAULT_COLOR_FORMAT = ColorFormat(sRGB)
STARTCOLORINDEX = 0
NCOLORS = 49
INPUT = "../run/peaks_data.csv"
OUTPUT = "../run/peaks_data.pdf"
PALETTE = "../data/yl_rd.ctb"

with open(PALETTE) as file:
    pal_raw = file.read()
pal = StackPalette()
for row in pal_raw.split("\n")[STARTCOLORINDEX:NCOLORS + 1]:
    if row:
        color = tuple(int(x) for x in row.split(" ")[1:])
        pal.add(color)

colorscale = []
for i, j in zip(range(0, len(pal) - 1), range(1, len(pal))):
    colorscale.append((i / (len(pal) - 1), pal[i].hex()))
    colorscale.append((j / (len(pal) - 1), pal[j].hex()))

df = pd.read_csv(INPUT, header=None)
max_food = df.values.max()
fig = make_subplots(len(df), 1)
for i in range(len(df)):
    row = df.iloc[i]
    trace = go.Scatter(
        x=list(range(1, len(row) + 1)),
        y=row,
        mode="markers",
        marker=go.scatter.Marker(color=row, colorscale=colorscale, cmax=max_food)
    )
    fig.add_trace(trace, row=i + 1, col=1)
fig.update_layout(showlegend=False, width=1000, height=1000)

fig.write_image(file=OUTPUT)

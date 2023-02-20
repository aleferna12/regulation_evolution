import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from plotly.subplots import make_subplots

POP = 500000
JMED = 14
JALPHA = 7
MULT = [1, 2, 3, 4, 5, 6]
MULT = np.array(MULT)


def main():
    jkeys = make_bitstrings(POP, len(MULT))
    jlocks = make_bitstrings(POP, len(MULT))
    ccJs = []
    gammas = []
    it = enumerate(zip(jkeys, jlocks))
    for i, (jkey1, jlock1) in it:
        for j, (jkey2, jlock2) in it:
            if i != j:
                ccJ = cell_cell_J(jkey1, jlock1, jkey2, jlock2)
                ccJs.append(ccJ)
                gammas.append(gamma(JMED, JALPHA, ccJ))

    fig = make_subplots(2, 1)
    values, counts = np.unique(ccJs, return_counts=True)
    fig.add_trace(go.Scatter(x=values, y=counts / np.sum(counts)), row=1, col=1)
    values, counts = np.unique(gammas, return_counts=True)
    fig.add_trace(go.Scatter(x=values, y=counts / np.sum(counts)), row=2, col=1)
    fig.update_layout(xaxis_title="Jcc",
                      yaxis_title="relative freq.",
                      title="Weights: " + " ".join(str(x) for x in MULT))
    fig.update_yaxes(range=[0, 0.25])
    fig.update_xaxes(range=[0, 50], row=1, col=1)
    fig.update_xaxes(range=[-25, 25], row=2, col=1)
    # fig.write_html(Path("~/Desktop/Jcc.html").expanduser())
    fig.show()


def make_bitstrings(n, length):
    rng = np.random.default_rng()
    return rng.integers(2, size=(n, length))


def cell_cell_J(k1, l1, k2, l2):
    j1 = np.sum((k1 == l2) * MULT)
    j2 = np.sum((k2 == l1) * MULT)
    return j1 + j2


def gamma(Jmed, Jalpha, Jcc):
    return Jmed - (Jalpha + Jcc) / 2


if __name__ == "__main__":
    main()
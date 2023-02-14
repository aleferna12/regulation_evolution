import plotly.express as px
import numpy as np

LEN = 12
POP = 500000


def main():
    jkeys = make_bitstrings(POP, LEN)
    jlocks = make_bitstrings(POP, LEN)
    ccJs = []
    it = enumerate(zip(jkeys, jlocks))
    for i, (jkey1, jlock1) in it:
        for j, (jkey2, jlock2) in it:
            if i != j:
                ccJs.append(cell_cell_J(jkey1, jlock1, jkey2, jlock2))

    values, counts = np.unique(ccJs, return_counts=True)
    fig = px.line(x=values, y=counts / np.sum(counts))
    fig.update_layout(xaxis_title="Jcc", yaxis_title="relative freq.")
    # fig.write_html(Path("~/Desktop/Jcc.html").expanduser())
    fig.show()


def make_bitstrings(n, length):
    rng = np.random.default_rng()
    return rng.integers(2, size=(n, length))


def cell_cell_J(k1, l1, k2, l2):
    mult = np.array([1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3])
    j1 = np.sum((k1 == l2) * mult)
    j2 = np.sum((k2 == l1) * mult)
    return j1 + j2


if __name__ == "__main__":
    main()
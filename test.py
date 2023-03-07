import numpy as np
import plotly.express as px
    

def fill_diag(a, fill, k=0):
    rowidx, colidx = np.arange(a.shape[1]), np.arange(a.shape[1])
    colidx = colidx.copy()  # rowidx and colidx share the same buffer
    if k > 0:
        colidx += k
        stop = min(a.shape[0], a.shape[1] - k)
    else:
        rowidx -= k
        stop = min(a.shape[1], a.shape[0] + k)
    a[rowidx[:stop], colidx[:stop]] = fill


img = np.zeros((100, 50))
for i in range(-100, 50, 5):
    fill_diag(img, 1, i)

px.imshow(img).show()
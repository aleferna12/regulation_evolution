"""This script distributes the area of a lattice to points that represent peaks of resources.

Inside the area around each peak is where the strength of the gradient is determined by that peak.

It's not very fast, but can be improved by using ray-tracing (this is done).

The other thing that can be done to improve is add points one by one and for each new point:
    1. Query which is the closest old point from the new point
    2. Add the closest point and new point to a list
    4. Ask what points do the closest point share classifiers with (interface list) and
    add those to the previously created list with new point and closest point
    5. For the new point, the closest point, and the points that bordered the closest point (all
    stored in our list):
        1. Re-calculate the polygon but only consider classifiers between the points in the list
        2. Recreate the interface list with the points that share classifiers with this point
"""

import plotly.graph_objects as go
import numpy as np
from math import dist, tan, pi, atan2
from scipy.spatial import ConvexHull
from colorir import *


def main():
    latt = (1000, 1000)
    n = 100
    # Number of rays used in ray-tracing
    # -1 disables it
    # Increasing this too much is worse than disabling it
    ray_n = 10
    # Float arithmetic needs to be corrected
    # A low precision may cause areas to overextend, while a high one may cause some areas not
    # to be assigned to any point
    precision = 4

    np.random.seed(2)
    # Centers of the grads
    points = np.random.rand(n, 2)
    points[:, 0] *= latt[0]
    points[:, 1] *= latt[1]
    points = [Point(*p) for p in points]

    fig = go.Figure()
    colors = Grad(
        StackPalette.load("spectral"),
        color_format=ColorFormat(
            sRGB,
            include_a=True
        )
    ).n_colors(n)

    ps = []
    for p in points:
        ps = add_to_latt(p, ps, latt, ray_n, precision)
    for i, point in enumerate(points):
        poly = get_poly(point, points, latt, ray_n, precision)
        poly = np.array(poly)
        # Theoretically this should always work but it depends on precision
        try:
            hull = ConvexHull(poly)
        except:
            hull = None

        color = colors[i]
        color = (color.r, color.g, color.b, 0.75)
        color = "rgba" + str(color)
        fig.add_trace(go.Scatter(
            x=[point[0]],
            y=[point[1]],
            marker_color=color
        ))
        if hull is not None:
            fig.add_trace(go.Scatter(
                mode="lines",
                x=poly[hull.vertices, 0],
                y=poly[hull.vertices, 1],
                line_color=color,
                fill="toself"
            ))
    fig.update_xaxes(range=[0, latt[0]], constrain="domain")
    fig.update_yaxes(range=[0, latt[1]], scaleanchor="x", scaleratio=1)
    fig.update_layout(showlegend=False)
    fig.show()


def add_to_latt(p, ps, latt, ray_n, prec):
    pc = closest_p(p, ps)
    if pc is not None:
        for pu in pc.infs:
            update_vertices(pu, [p, pc] + pu.infs, latt, prec)
        update_vertices(p, [pc] + pc.infs, latt, prec)
        update_vertices(pc, [p] + pc.infs, latt, prec)
    return ps + [p]


def update_vertices(p, check_ps, latt, prec):
    # p.infs, ncs = neighbours(p, check_ps, prec)
    # p.vs = find_all_vertices(p, ncs, latt, prec)
    pass


def get_poly(p, ps, latt, ray_n, prec):
    """Check all classifiers. Slow but precise."""
    cs = []
    for p2 in ps:
        c = get_classifier(p, p2)
        if c is not None:
            cs.append(c)

    if ray_n > 0:
        ray_cs = []
        d_step = 2 * pi / ray_n
        for i in range(ray_n):
            angle = i * d_step
            m = tan(angle)
            n = p[1] - p[0] * m
            fc = ray_trace(p, (m, n), angle, cs)
            if fc is not None:
                ray_cs.append(fc)
        cs = ray_cs

    return [v[1] for v in get_all_vertices(p, ps, cs, latt, prec)]


def closest_p(p, ps):
    closest = (None, float("inf"))
    for p2 in ps:
        d = dist(p, p2)
        if d < closest[1]:
            closest = (p2, d)
    return closest[0]


def ray_trace(p, q, angle, cs):
    """Selects the closest classifier in cs that intersects q"""
    closest = (None, float("inf"))
    q_quad = get_quad(angle)
    for c in cs:
        p2 = get_intersection(q, c)

        d = dist(p, p2)
        if d < closest[1]:
            x = p[0] - p2[0]
            y = p[1] - p2[1]
            angle2 = atan2(y, x)
            if get_quad(angle2) == q_quad:
                closest = (c, d)
    return closest[0]


def get_quad(angle):
    angle = angle % (2 * pi)
    if 0 <= angle < pi / 2:
        return 0
    if pi / 2 <= angle < pi:
        return 1
    if pi <= angle < 3 * pi / 2:
        return 2
    return 3


def get_classifier(p1, p2):
    if p1 == p2:
        return None

    r = get_line(p1, p2)
    # Average point
    p3 = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
    return get_perpendicular(p3, r[0])


def get_all_vertices(p, ps, cs, latt, prec):
    # For every classifier we need to query whether it has "direct line of sight" to p
    # Vertices are the putative points that match this condition
    # They will later be filtered
    vertices = [(0, 0), (0, latt[1]), (latt[0], 0), (latt[0], latt[1])]
    vertices = list(zip([None] * len(vertices), vertices))
    for p2, c in zip(ps, cs):
        if c is None:
            continue
        vertices += list(zip([None] * 4, get_intersection_latt(c, latt)))

        for o_c in cs:
            if o_c == c or o_c is None:
                continue
            pv = get_intersection(c, o_c)
            if in_bounds(pv, latt):
                vertices.append((p2, pv))
    return [v for v in vertices if not split_by_classif(p, v[1], cs, prec)]


def get_box(poly):
    return [
        min(pp[0] for pp in poly),
        min(pp[1] for pp in poly),
        max(pp[0] for pp in poly),
        max(pp[1] for pp in poly)
    ]


def validate(vs, lower, upper, repl):
    ret = [v if v is not None and lower < v < upper else repl for v in vs]
    return ret


def in_bounds(p, latt):
    return 0 <= p[0] <= latt[0] and 0 <= p[1] <= latt[1]


def split_by_classif(p1, p2, cs, prec):
    """Whether the line from p1 to p2 is split by any classifier in cs."""
    # Rounding is required for float error not to take place
    r = get_line(p1, p2)
    p1 = (round(p1[0], prec), round(p1[1], prec))
    p2 = (round(p2[0], prec), round(p2[1], prec))
    for c in cs:
        pc = get_intersection(r, c)
        pc = (round(pc[0], prec), round(pc[1], prec))
        if min(p1[0], p2[0]) < pc[0] < max(p1[0], p2[0]) and min(p1[1], p2[1]) < pc[1] < max(p1[1], p2[1]):
            return True
    return False


def get_intersection_latt(c, latt):
    return [(0, c[1]),
            (-c[1] / c[0], 0),
            (latt[0], latt[0] * c[0] + c[1]),
            ((latt[1] - c[1]) / c[0], latt[1])]


def get_intersection(c1, c2):
    if c1[0] == c2[0]:
        return None
    x = (c1[1] - c2[1]) / (c2[0] - c1[0])
    y = (c1[1] * c2[0] - c2[1] * c1[0]) / (c2[0] - c1[0])
    return x, y


def box_from_coords(ps):
    # minx, miny, maxx, maxy
    return (
        ps[0][0],
        ps[1][1],
        ps[2][0],
        ps[3][1]
    )


# This is similar but not the same as get_line
def get_perpendicular(pm, m):
    mp = -1 / m
    np_ = pm[1] - pm[0] * mp
    return mp, np_


def get_line(p1, p2):
    m = (p2[1] - p1[1]) / (p2[0] - p1[0])
    n = (p1[0] * p2[1] - p2[0] * p1[1]) / (p1[0] - p2[0])
    return m, n


def lattice_ref_points(p, latt):
    prs = [
        (0, p[1]),
        (p[0], 0),
        (latt[0], p[1]),
        (p[0], latt[1])
    ]
    return [dist(p, pr) for pr in prs]


class Point:
    def __init__(self, x, y, infs=None, vs=None):
        self.x, self.y = x, y
        self.infs = [] if infs is None else infs
        self.vs = [] if infs is None else vs

    def __getitem__(self, item):
        if item == 0:
            return self.x
        elif item == 1:
            return self.y
        raise Exception("2d point has only 2 dimensions (duhh)")

    def __iter__(self):
        return iter((self.x, self.y))

    def __repr__(self):
        return f"{self.x}, {self.y}"


if __name__ == "__main__":
    main()

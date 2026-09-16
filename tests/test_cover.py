import numpy as np
from enclosing.geometry import make_polygon, lower_bounds
from enclosing.refine import iterative_cover


def test_square_small():
    poly = make_polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    R = 1.0
    res = iterative_cover(poly, R, h0=0.5, s0=0.5, max_iter=6,
                          solver="auto", verbose=False)
    assert res is not None
    assert res["check"]["covered"], res["history"]
    # cota de Kershner: ceil(4 / 2.598) = 2
    assert lower_bounds(poly, R)["kershner_bound"] == 2
    assert 2 <= res["k"] <= 8


def test_l_shape():
    poly = make_polygon([(0, 0), (3, 0), (3, 1), (1, 1), (1, 3), (0, 3)])
    R = 0.8
    res = iterative_cover(poly, R, h0=0.4, s0=0.4, max_iter=6,
                          solver="auto", verbose=False)
    assert res is not None
    assert res["check"]["covered"], res["history"]
    assert res["k"] >= res["bounds"]["kershner_bound"]

import numpy as np
from enclosing.geometry import make_polygon, lower_bounds, shape_p
from enclosing.sampling import candidate_centers
from enclosing.cover import build_coverage
from enclosing.verify import check_coverage
from enclosing.refine import iterative_cover


def test_shape_p():
    assert shape_p("circle") == 2
    assert shape_p("square") == np.inf
    try:
        shape_p("triangle")
    except ValueError:
        pass
    else:
        raise AssertionError("shape invalido no lanza ValueError")


def test_square_check_exact():
    # cuatro cuadrados s=1 en (1,1),(1,3),(3,1),(3,3) cubren el 4x4 exacto
    poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
    C = np.array([[1.0, 1.0], [1.0, 3.0], [3.0, 1.0], [3.0, 3.0]])
    chk = check_coverage(poly, C, 1.0, shape="square")
    assert chk["covered"]
    assert chk["uncovered_ratio"] == 0.0
    assert chk["uncovered_area"] == 0.0
    # con discos R=1 los mismos centros dejan huecos
    chk_c = check_coverage(poly, C, 1.0, shape="circle")
    assert not chk_c["covered"]


def test_square_coverage_chebyshev():
    # norma infinito: (0.9, 0.9) esta a <=1 del origen; (1.1, 0) no
    dem = np.array([[0.9, 0.9], [1.1, 0.0]])
    cand = np.array([[0.0, 0.0]])
    cl, ui = build_coverage(dem, cand, 1.0, shape="square")
    assert cl[0] == [0]
    assert ui == [1]
    # en euclidea (0.9,0.9) esta a 1.27 > 1: nadie la cubre
    cl2, ui2 = build_coverage(dem, cand, 1.0, shape="circle")
    assert ui2 == [0, 1]


def test_square_candidates_prune():
    poly = make_polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    cand = candidate_centers(poly, 1.0, 1.0, method="square", shape="square")
    assert len(cand) > 0
    # todo candidato toca al poligono con su cuadrado
    from shapely.geometry import box
    for x, y in cand:
        assert box(x - 1.0, y - 1.0, x + 1.0, y + 1.0).intersects(poly)


def test_square_small():
    poly = make_polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    s = 1.0
    b = lower_bounds(poly, s, shape="square")
    assert b["area_bound"] == 1  # ceil(4/4)
    assert b["kershner_bound"] == 1
    res = iterative_cover(poly, s, h0=0.5, s0=0.5, max_iter=6,
                          solver="auto", verbose=False, shape="square")
    assert res is not None
    assert res["check"]["covered"], res["history"]
    assert res["check"]["uncovered_ratio"] == 0.0
    assert 1 <= res["k"] <= 6


def test_square_l_shape():
    poly = make_polygon([(0, 0), (3, 0), (3, 1), (1, 1), (1, 3), (0, 3)])
    s = 0.8
    res = iterative_cover(poly, s, h0=0.4, s0=0.4, max_iter=6,
                          solver="auto", verbose=False, shape="square")
    assert res is not None
    assert res["check"]["covered"], res["history"]
    assert res["k"] >= res["bounds"]["area_bound"]


def test_structural_seeds():
    # el vertice (0,0) genera, entre otras, la semilla (-1,-1) y (1,1);
    # (-1,-1) toca al poligono con su cuadrado y debe conservarse
    poly = make_polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
    cand = candidate_centers(poly, 1.0, 1.0, method="square",
                             shape="square", seed="structural")
    assert any(np.allclose(c, [-1.0, -1.0]) for c in cand)
    grid = candidate_centers(poly, 1.0, 1.0, method="square",
                             shape="square", seed="grid")
    assert len(cand) >= len(grid)
    # con semillas estructurales el 2x2 s=1 se cubre con k pequeno y exacto
    res = iterative_cover(poly, 1.0, h0=0.5, s0=1.0, max_iter=4,
                          solver="auto", verbose=False, shape="square",
                          seed="structural")
    assert res["check"]["covered"]
    assert res["check"]["uncovered_ratio"] == 0.0
    assert res["k"] <= 4

"""Matriz de cobertura demanda x candidatos."""
import numpy as np
from scipy.spatial import cKDTree
from .geometry import shape_p


def build_coverage(demand, candidates, R, tol=1e-9, shape="circle"):
    """Devuelve cover_lists[i] = lista de candidatos que cubren demanda i.

    shape="circle": |p-c|_2 <= R. shape="square" (R = semilado, ejes
    paralelos): |p-c|_inf <= R. El MILP binario no cambia: solo la matriz.

    Tambien indices de demandas no cubiertas por ningun candidato
    (infactibilidad del discretizado).
    """
    p = shape_p(shape)
    Nd = len(demand)
    Nc = len(candidates)
    if Nc == 0:
        return [[] for _ in range(Nd)], list(range(Nd))
    tree = cKDTree(candidates)
    ball = tree.query_ball_point(demand, r=R + tol, p=p)
    cover_lists = [sorted(map(int, lst)) for lst in ball]
    uncovered = [i for i, lst in enumerate(cover_lists) if len(lst) == 0]
    return cover_lists, uncovered

"""Enclosing: cobertura de regiones planas con circulos iguales (caso 1)."""
from .geometry import make_polygon, lower_bounds
from .sampling import demand_points, candidate_centers
from .cover import build_coverage
from .solvers import solve_set_cover
from .verify import check_coverage
from .refine import iterative_cover

__all__ = [
    "make_polygon",
    "lower_bounds",
    "demand_points",
    "candidate_centers",
    "build_coverage",
    "solve_set_cover",
    "check_coverage",
    "iterative_cover",
]

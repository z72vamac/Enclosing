"""Utilidades geometricas basicas."""
import math
import numpy as np
from shapely.geometry import Polygon


def make_polygon(coords, holes=None):
    """Crea un Polygon valido a partir de coordenadas."""
    poly = Polygon(coords, holes=holes)
    if not poly.is_valid:
        poly = poly.buffer(0)
    if poly.is_empty or not poly.is_valid:
        raise ValueError("Poligono invalido incluso tras buffer(0).")
    return poly


def hex_cell_area(R):
    """Area del hexagono regular inscrito en circulo de radio R."""
    return 3.0 * math.sqrt(3.0) / 2.0 * R * R


def lower_bounds(poly, R):
    """Cotas inferiores del numero de circulos.

    area_bound: ceil(A / pi R^2) (sin solape, inalcanzable).
    kershner: ceil(A / A_hex) (asintotico hexagonal, mas fuerte).
    """
    A = float(poly.area)
    area_b = math.ceil(A / (math.pi * R * R) - 1e-12)
    kersh = math.ceil(A / hex_cell_area(R) - 1e-12)
    return {"area": A, "perimeter": float(poly.length),
            "area_bound": max(area_b, 1), "kershner_bound": max(kersh, 1)}


def bounds_expanded(poly, delta):
    minx, miny, maxx, maxy = poly.bounds
    return (minx - delta, miny - delta, maxx + delta, maxy + delta)

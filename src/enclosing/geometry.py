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


def shape_p(shape):
    """Orden de Minkowski asociado a la forma: 2 (circulo) o inf (cuadrado).

    Con shape="square", R se interpreta como semilado s (orientacion fija).
    """
    if shape == "circle":
        return 2
    if shape == "square":
        return np.inf
    raise ValueError(f"shape desconocido: {shape!r} (usa 'circle' o 'square')")


def lower_bounds(poly, R, shape="circle"):
    """Cotas inferiores del numero de piezas.

    circle: area_bound = ceil(A / pi R^2); kershner con hexagono inscrito.
    square (semilado R): area_bound = ceil(A / 4R^2); los cuadrados teselan
    el plano, asi que la cota de area es el analogo de Kershner.
    """
    A = float(poly.area)
    if shape == "square":
        tile = 4.0 * R * R
        area_b = math.ceil(A / tile - 1e-12)
        return {"area": A, "perimeter": float(poly.length),
                "area_bound": max(area_b, 1), "kershner_bound": max(area_b, 1)}
    if shape != "circle":
        raise ValueError(f"shape desconocido: {shape!r}")
    area_b = math.ceil(A / (math.pi * R * R) - 1e-12)
    kersh = math.ceil(A / hex_cell_area(R) - 1e-12)
    return {"area": A, "perimeter": float(poly.length),
            "area_bound": max(area_b, 1), "kershner_bound": max(kersh, 1)}


def bounds_expanded(poly, delta):
    minx, miny, maxx, maxy = poly.bounds
    return (minx - delta, miny - delta, maxx + delta, maxy + delta)

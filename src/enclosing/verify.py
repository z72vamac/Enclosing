"""Verificacion continua con shapely: union de piezas vs poligono."""
from shapely.geometry import Point, box
from shapely.ops import unary_union


def check_coverage(poly, centers, R, buffer_resolution=64, shape="circle"):
    """Comprueba si la union de piezas cubre poly.

    shape="circle": discos de radio R (buffer poligonal aproximado).
    shape="square": cuadrados ejes-paralelos de semilado R (exactos).

    Devuelve dict con covered, uncovered_area, uncovered_ratio, uncovered_geom.
    """
    if shape not in ("circle", "square"):
        raise ValueError(f"shape desconocido: {shape!r}")
    if len(centers) == 0:
        return {"covered": False, "uncovered_area": float(poly.area),
                "uncovered_ratio": 1.0, "uncovered_geom": poly,
                "union": None}
    if shape == "square":
        pieces = [box(float(x) - float(R), float(y) - float(R),
                      float(x) + float(R), float(y) + float(R))
                  for x, y in centers]
    else:
        pieces = [Point(float(x), float(y)).buffer(float(R),
                  resolution=int(buffer_resolution)) for x, y in centers]
    union = unary_union(pieces)
    uncovered = poly.difference(union)
    ua = float(uncovered.area) if not uncovered.is_empty else 0.0
    ratio = ua / float(poly.area) if poly.area > 0 else 0.0
    # circulo: el buffer es una aproximacion poligonal del disco,
    # deja rendijas de ~R^2/resolution^2. 1e-7 es conservador.
    # cuadrado: geometria exacta; queda solo redondeo flotante.
    tol = 1e-9 if shape == "square" else 1e-7
    covered = bool(uncovered.is_empty or ratio < tol)
    return {"covered": covered, "uncovered_area": ua,
            "uncovered_ratio": ratio, "uncovered_geom": uncovered,
            "union": union}

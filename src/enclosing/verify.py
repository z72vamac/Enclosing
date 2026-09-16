"""Verificacion continua con shapely: union de discos vs poligono."""
from shapely.geometry import Point
from shapely.ops import unary_union


def check_coverage(poly, centers, R, buffer_resolution=64):
    """Comprueba si la union de discos cubre poly.

    Devuelve dict con covered, uncovered_area, uncovered_ratio, uncovered_geom.
    """
    if len(centers) == 0:
        return {"covered": False, "uncovered_area": float(poly.area),
                "uncovered_ratio": 1.0, "uncovered_geom": poly,
                "union": None}
    disks = [Point(float(x), float(y)).buffer(float(R),
             resolution=int(buffer_resolution)) for x, y in centers]
    union = unary_union(disks)
    uncovered = poly.difference(union)
    ua = float(uncovered.area) if not uncovered.is_empty else 0.0
    ratio = ua / float(poly.area) if poly.area > 0 else 0.0
    # tolerancia numerica: el buffer es una aproximacion poligonal del disco,
    # deja rendijas de ~R^2/resolution^2. 1e-7 es conservador.
    covered = bool(uncovered.is_empty or ratio < 1e-7)
    return {"covered": covered, "uncovered_area": ua,
            "uncovered_ratio": ratio, "uncovered_geom": uncovered,
            "union": union}

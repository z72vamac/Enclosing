"""Muestreo de puntos de demanda y centros candidatos."""
import numpy as np
from shapely.geometry import Point
from shapely.prepared import prep
from .geometry import bounds_expanded


def _filter_inside(points, poly, tol=1e-9):
    """Filtra puntos dentro de poly (incluye borde)."""
    if len(points) == 0:
        return points
    pre = prep(poly)
    # buffer tiny para incluir borde con tolerancia
    buf = poly.buffer(tol)
    pre_buf = prep(buf)
    mask = np.array([pre_buf.covers(Point(p)) for p in points], dtype=bool)
    return points[mask]


def _square_lattice(x0, y0, x1, y1, s):
    xs = np.arange(x0, x1 + 0.5 * s, s)
    ys = np.arange(y0, y1 + 0.5 * s, s)
    xx, yy = np.meshgrid(xs, ys)
    return np.column_stack([xx.ravel(), yy.ravel()])


def _hex_lattice(x0, y0, x1, y1, s):
    """Red hexagonal: distancia entre centros vecinos = s."""
    dy = s * np.sqrt(3.0) / 2.0
    ys = np.arange(y0, y1 + 0.5 * s, dy)
    pts = []
    for j, y in enumerate(ys):
        offset = (s / 2.0) if (j % 2 == 1) else 0.0
        xs = np.arange(x0 + offset, x1 + 0.5 * s, s)
        for x in xs:
            pts.append((x, y))
    if not pts:
        return np.zeros((0, 2))
    return np.array(pts, dtype=float)


def interior_points(poly, spacing, method="square"):
    minx, miny, maxx, maxy = poly.bounds
    if method == "hex":
        grid = _hex_lattice(minx, miny, maxx, maxy, spacing)
    else:
        grid = _square_lattice(minx, miny, maxx, maxy, spacing)
    return _filter_inside(grid, poly)


def boundary_points(poly, spacing):
    """Muestrea exterior e interiores cada `spacing`."""
    pts = []
    for ring in [poly.exterior, *poly.interiors]:
        L = ring.length
        n = max(int(round(L / spacing)), 1)
        for i in range(n):
            p = ring.interpolate(i / n * L)
            pts.append((p.x, p.y))
    if not pts:
        return np.zeros((0, 2))
    return np.array(pts, dtype=float)


def demand_points(poly, h, hb=None, method="square"):
    """Demanda = interior (paso h) + borde (paso hb=h/2 por defecto)."""
    if hb is None:
        hb = h / 2.0
    inner = interior_points(poly, h, method=method)
    bord = boundary_points(poly, hb)
    if len(inner) == 0:
        allp = bord
    elif len(bord) == 0:
        allp = inner
    else:
        allp = np.vstack([inner, bord])
    # deduplicar por redondeo
    if len(allp) == 0:
        # poligono diminuto: usa un punto representativo
        p = poly.representative_point()
        return np.array([[p.x, p.y]])
    key = np.round(allp / (hb / 4.0)).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    return allp[np.sort(idx)]


def candidate_centers(poly, R, spacing, method="hex"):
    """Candidatos en bbox expandido por R, podados a dist(poly) <= R."""
    x0, y0, x1, y1 = bounds_expanded(poly, R)
    if method == "hex":
        grid = _hex_lattice(x0, y0, x1, y1, spacing)
    else:
        grid = _square_lattice(x0, y0, x1, y1, spacing)
    if len(grid) == 0:
        return grid
    # poda: descarta candidatos que no tocan P expandido por R
    buf = poly.buffer(R + 1e-9)
    pre = prep(buf)
    mask = np.array([pre.covers(Point(p)) for p in grid], dtype=bool)
    return grid[mask]


def sample_points_in_geom(geom, spacing):
    """Muestrea puntos dentro de geom (Polygon o MultiPolygon)."""
    if geom.is_empty:
        return np.zeros((0, 2))
    geoms = [geom] if geom.geom_type == "Polygon" else list(geom.geoms)
    out = []
    for g in geoms:
        if g.is_empty or g.area <= 0:
            # agujero degenerado: punto representativo
            try:
                p = g.representative_point()
                out.append((p.x, p.y))
            except Exception:
                pass
            continue
        minx, miny, maxx, maxy = g.bounds
        grid = _square_lattice(minx, miny, maxx, maxy, spacing)
        pre = prep(g.buffer(1e-12))
        for p in grid:
            if pre.covers(Point(p)):
                out.append((p[0], p[1]))
        # garantiza al menos un punto por componente
        p = g.representative_point()
        out.append((p.x, p.y))
    if not out:
        return np.zeros((0, 2))
    arr = np.array(out, dtype=float)
    key = np.round(arr / (spacing / 4.0 + 1e-12)).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    return arr[np.sort(idx)]

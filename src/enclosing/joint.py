"""Variante A: recubrimiento + ruta en uno (joint location-routing).

Dos piezas:
1. solve_joint_exact: MILP conjunto (cobertura + TSP/MTZ) sobre candidatos
   discretos. Solo instancias pequenas (m ~ <=60).
2. alternating_refine: heuristica para instancias grandes. Alterna
   re-posicionar centros (con cobertura de demanda fija) y re-optimizar tour.
"""
import numpy as np
from scipy.spatial import cKDTree
from scipy.optimize import minimize
from .geometry import shape_p
from .tsp import dist_matrix, tour_length, two_opt, solve_tsp_heuristic


# ---------------- exacto ----------------

def _extract_tour_from_x(xval, selected, dep):
    succ = {}
    for i in selected:
        for j in selected:
            if xval.get((i, j), 0) > 0.5:
                succ[i] = j
                break
    tour = [dep]
    for _ in range(len(selected) - 1):
        nxt = succ.get(tour[-1])
        if nxt is None:
            return None
        tour.append(nxt)
    return tour if succ.get(tour[-1]) == dep else None


def _find_subtours(xval, selected):
    succ = {}
    for i in selected:
        for j in selected:
            if xval.get((i, j), 0) > 0.5:
                succ[i] = j
                break
    seen, cycles = set(), []
    for s in selected:
        if s in seen:
            continue
        cyc, cur = [], s
        while cur not in seen and cur in succ:
            seen.add(cur)
            cyc.append(cur)
            cur = succ[cur]
        if cyc and cur == s:
            cycles.append(cyc)
    return cycles

def solve_joint_dfj(candidates, cover_lists, depot_idx=0, time_limit=300,
                    warm_y=None, warm_tour=None, output_flag=0):
    """MILP conjunto con DFJ perezoso: min tour s.a. cobertura + TSP.

    Master: cobertura + grados ligados. Los subtours se eliminan con
    lazy cuts en incumbentes enteros. Escala mucho mejor que MTZ.
    """
    import gurobipy as gp
    from gurobipy import GRB
    C = np.asarray(candidates, dtype=float)
    m = len(C)
    D = dist_matrix(C)
    dep = int(depot_idx)
    mod = gp.Model()
    mod.setParam("OutputFlag", output_flag)
    mod.setParam("LazyConstraints", 1)
    if time_limit is not None:
        mod.setParam("TimeLimit", float(time_limit))
    y = mod.addVars(m, vtype="B", name="y")
    x = mod.addVars(m, m, vtype="B", name="x")
    for i in range(m):
        x[i, i].UB = 0
    for i, lst in enumerate(cover_lists):
        if len(lst) == 0:
            return {"feasible": False, "message": f"demanda {i} sin candidatos"}
        mod.addConstr(gp.quicksum(y[j] for j in lst) >= 1)
    for i in range(m):
        mod.addConstr(gp.quicksum(x[i, j] for j in range(m)) == y[i])
        mod.addConstr(gp.quicksum(x[j, i] for j in range(m)) == y[i])
    mod.addConstr(y[dep] == 1)
    mod.setObjective(gp.quicksum(D[i, j] * x[i, j]
                                 for i in range(m) for j in range(m)))
    if warm_y is not None and warm_tour is not None:
        for j in warm_y:
            y[j].Start = 1.0
        t = list(warm_tour) + [warm_tour[0]]
        for a, b in zip(t[:-1], t[1:]):
            x[a, b].Start = 1.0

    def cb(model, where):
        if where != GRB.Callback.MIPSOL:
            return
        yv = model.cbGetSolution(y)
        xv = model.cbGetSolution(x)
        sel = [i for i in range(m) if yv[i] > 0.5]
        if len(sel) <= 2:
            return
        xval = {(i, j): xv[i, j] for i in sel for j in sel}
        for cyc in _find_subtours(xval, sel):
            if len(cyc) < len(sel):
                model.cbLazy(gp.quicksum(x[i, j] for i in cyc for j in cyc)
                             <= len(cyc) - 1)

    mod.optimize(cb)
    if mod.SolCount == 0:
        return {"feasible": False, "status": mod.status}
    yy = sorted(j for j in range(m) if y[j].X > 0.5)
    xval = {(i, j): x[i, j].X for i in yy for j in yy}
    tour = _extract_tour_from_x(xval, yy, dep)
    if tour is None:  # incumbent con subtours al agotar tiempo
        cycs = _find_subtours(xval, yy)
        tour = max(cycs, key=len) if cycs else yy
    L = tour_length(C, tour, closed=True)
    gap = getattr(mod, "MIPGap", None)
    try:
        gap = float(gap)
    except Exception:
        gap = None
    return {"feasible": True, "selected": yy, "tour": tour,
            "length": L, "k": len(yy), "status": mod.status,
            "gap": gap,
            "optimal": bool(mod.status == 2 and (gap or 0) < 1e-6)}


def solve_joint_exact(*args, **kwargs):
    """Alias de compatibilidad: redirige al DFJ (escala mejor que MTZ)."""
    return solve_joint_dfj(*args, **kwargs)


def solve_joint_mtz(candidates, cover_lists, depot_idx=0, time_limit=300,
                    warm_y=None, warm_tour=None, output_flag=0):
    """Formulacion MTZ clasica (referencia: relaja peor que DFJ)."""
    import gurobipy as gp
    C = np.asarray(candidates, dtype=float)
    m = len(C)
    D = dist_matrix(C)
    dep = int(depot_idx)
    mod = gp.Model()
    mod.setParam("OutputFlag", output_flag)
    if time_limit is not None:
        mod.setParam("TimeLimit", float(time_limit))
    y = mod.addVars(m, vtype="B", name="y")
    x = mod.addVars(m, m, vtype="B", name="x")
    for i in range(m):
        x[i, i].UB = 0
    for i, lst in enumerate(cover_lists):
        if len(lst) == 0:
            return {"feasible": False, "message": f"demanda {i} sin candidatos"}
        mod.addConstr(gp.quicksum(y[j] for j in lst) >= 1)
    for i in range(m):
        mod.addConstr(gp.quicksum(x[i, j] for j in range(m)) == y[i])
        mod.addConstr(gp.quicksum(x[j, i] for j in range(m)) == y[i])
    u = mod.addVars(m, vtype="C", lb=0, ub=m - 1, name="u")
    mod.addConstr(y[dep] == 1)
    mod.addConstr(u[dep] == 0)
    for i in range(m):
        if i == dep:
            continue
        for j in range(m):
            if j == dep or j == i:
                continue
            mod.addConstr(u[i] - u[j] + m * x[i, j] <= m - 1)
    mod.setObjective(gp.quicksum(D[i, j] * x[i, j]
                                 for i in range(m) for j in range(m)))
    if warm_y is not None and warm_tour is not None:
        for j in warm_y:
            y[j].Start = 1.0
        t = list(warm_tour) + [warm_tour[0]]
        for a, b in zip(t[:-1], t[1:]):
            x[a, b].Start = 1.0
        for k, node in enumerate(warm_tour):
            if node != dep:
                u[node].Start = float(k)
    mod.optimize()
    if mod.SolCount == 0:
        return {"feasible": False, "status": mod.status}
    yy = sorted(j for j in range(m) if y[j].X > 0.5)
    xval = {(i, j): x[i, j].X for i in yy for j in yy}
    tour = _extract_tour_from_x(xval, yy, dep)
    if tour is None:
        cycs = _find_subtours(xval, yy)
        tour = max(cycs, key=len) if cycs else yy
    L = tour_length(C, tour, closed=True)
    try:
        gap = float(mod.MIPGap)
    except Exception:
        gap = None
    return {"feasible": True, "selected": yy, "tour": tour,
            "length": L, "k": len(yy), "status": mod.status,
            "gap": gap,
            "optimal": bool(mod.status == 2 and (gap or 0) < 1e-6)}


# ---------------- alternante ----------------

def _assign_demands(demands, centers, R, shape="circle"):
    """Cada demanda al centro mas cercano de entre los que la cubren."""
    dem = np.asarray(demands, dtype=float)
    C = np.asarray(centers, dtype=float)
    tree = cKDTree(C)
    asg = [[] for _ in range(len(C))]
    for p in dem:
        idxs = tree.query_ball_point(p, r=R + 1e-9, p=shape_p(shape))
        if not idxs:
            return None  # demanda descubierta
        j = min(idxs, key=lambda j: float(np.linalg.norm(C[j] - p)))
        asg[j].append(p)
    return [np.array(a) if len(a) else np.zeros((0, 2)) for a in asg]


def _demands_covered(demands, centers, R, shape="circle"):
    if len(centers) == 0:
        return len(demands) == 0
    tree = cKDTree(np.asarray(centers, dtype=float))
    n = tree.query_ball_point(np.asarray(demands, dtype=float), r=R + 1e-9,
                              p=shape_p(shape))
    return all(len(v) > 0 for v in n)


def _best_position(a, b, pts, R, x0, max_step=None, shape="circle"):
    """min |p-a|+|p-b| s.a. |p-q|<=R (SLSQP). Si falla, devuelve x0.

    shape="square": la cobertura es |p-q|_inf <= R, o sea restricciones de
    caja lineales en p (orientacion fija, sin no linealidad geometrica).
    max_step: radio de confianza alrededor de x0 (preserva cobertura
    continua entre rondas; la verificacion por ronda decide).
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    P = np.asarray(pts, dtype=float)
    x0 = np.asarray(x0, dtype=float)
    ord_ = shape_p(shape)

    def fun(p):
        return float(np.linalg.norm(p - a) + np.linalg.norm(p - b))

    def jac(p):
        g = np.zeros(2)
        for z in (a, b):
            d = p - z
            n = np.linalg.norm(d)
            if n > 1e-12:
                g = g + d / n
        return g

    cons = [{"type": "ineq",
             "fun": (lambda q: (lambda p: float(R - np.linalg.norm(p - q, ord=ord_))))(q)}
            for q in P]
    if max_step is not None:
        cons.append({"type": "ineq",
                     "fun": (lambda p: float(max_step - np.linalg.norm(p - x0)))})
    try:
        res = minimize(fun, x0, jac=jac, method="SLSQP",
                       constraints=cons,
                       options={"maxiter": 200, "ftol": 1e-12})
    except Exception:
        return x0
    if (res.success and fun(res.x) < fun(x0) - 1e-9
            and all(float(np.linalg.norm(res.x - q, ord=ord_)) <= R + 1e-7
                    for q in P)):
        return res.x
    return x0


def alternating_refine(poly, R, centers, tour, demands, n_rounds=5,
                       try_drops=True, verbose=True, verify_tol=1e-6,
                       max_step=None, shape="circle"):
    """Fija orden -> mueve centros (cobertura preservada) -> re-opt tour.

    Tras cada ronda verifica cobertura CONTINUA; si se rompe, revierte a la
    mejor ronda cubierta. `demands` debe ser una malla DENSA de guarda
    (p.ej. h=R/4): los movimientos solo preservan lo que la demanda impone.
    max_step (defecto R/8): radio de confianza por ronda.
    shape: "circle" o "square" (semilado R, restricciones de caja lineales).
    """
    from .verify import check_coverage
    if max_step is None:
        max_step = R / 8.0
    C = [np.asarray(c, dtype=float) for c in centers]
    t = list(tour)
    dem = np.asarray(demands, dtype=float)
    L = tour_length(np.array(C), t, closed=True)
    hist = [{"round": -1, "k": len(C), "length": L}]
    best = {"centers": np.array(C), "tour": list(t), "length": L,
            "k": len(C)}
    if verbose:
        print(f"[alt 0] k={len(C)} L={L:.4f}")

    for rd in range(int(n_rounds)):
        # 1. drops: elimina centros redundantes sobre la demanda
        if try_drops:
            order = sorted(range(len(C)),
                           key=lambda v: float(np.linalg.norm(C[v] - C[t[0]])))
            for v in order:
                if len(C) <= 2:
                    break
                keep = [i for i in range(len(C)) if i != v]
                if _demands_covered(dem, [C[i] for i in keep], R,
                                    shape=shape):
                    # remapea tour
                    t = [i if i < v else i - 1 for i in t if i != v]
                    C.pop(v)
        # 2. re-asigna y mueve cada centro hacia su segmento del tour
        asg = _assign_demands(dem, np.array(C), R, shape=shape)
        if asg is None:
            if verbose:
                print(f"[alt {rd}] demanda descubierta, paro")
            break
        arr = np.array(C)
        for pos, node in enumerate(t):
            if len(asg[node]) == 0:
                continue
            prev = arr[t[pos - 1]]
            nxt = arr[t[(pos + 1) % len(t)]]
            arr[node] = _best_position(prev, nxt, asg[node], R, arr[node],
                                       max_step=max_step, shape=shape)
        C = [r for r in arr]
        # 3. re-optimiza el orden
        D = dist_matrix(np.array(C))
        t = two_opt(t, D)
        L = tour_length(np.array(C), t, closed=True)
        chk_rd = check_coverage(poly, np.array(C), R, shape=shape)
        ok = bool(chk_rd["covered"]
                  or chk_rd["uncovered_ratio"] < verify_tol)
        hist.append({"round": rd, "k": len(C), "length": L,
                     "covered": ok,
                     "uncovered_ratio": chk_rd["uncovered_ratio"]})
        if verbose:
            print(f"[alt {rd+1}] k={len(C)} L={L:.4f} "
                  f"cubierto={ok} ratio={chk_rd['uncovered_ratio']:.2e}")
        if ok:
            best = {"centers": np.array(C), "tour": list(t),
                    "length": L, "k": len(C)}
        else:
            if verbose:
                print(f"[alt {rd+1}] revierto: se rompio cobertura continua")
            C = [r for r in best["centers"]]
            t = list(best["tour"])
            break
    chk = check_coverage(poly, best["centers"], R, shape=shape)
    return {"centers": best["centers"], "tour": best["tour"],
            "length": best["length"], "k": best["k"],
            "history": hist, "check": chk}


# ---------------- Voronoi continuo (mueve centros, cobertura exacta) ---

def _finite_voronoi_regions(vor, radius):
    """Receta estandar de scipy: cierra regiones infinitas a `radius`."""
    import numpy as np
    new_regions = []
    new_vertices = vor.vertices.tolist()
    center = vor.points.mean(axis=0)
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]
        if all(v >= 0 for v in vertices):
            new_regions.append(vertices)
            continue
        ridges = all_ridges.get(p1, [])
        new_region = [v for v in vertices if v >= 0]
        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:  # arista finita: ya incluida
                continue
            tangent = vor.points[p2] - vor.points[p1]
            tangent = tangent / np.linalg.norm(tangent)
            normal = np.array([-tangent[1], tangent[0]])
            midpoint = (vor.points[p1] + vor.points[p2]) / 2.0
            direction = np.sign(np.dot(midpoint - center, normal)) * normal
            far_point = vor.vertices[v2] + direction * radius
            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:, 1] - c[1], vs[:, 0] - c[0])
        new_region = [new_region[i] for i in np.argsort(angles)]
        new_regions.append(new_region)
    return new_regions, np.asarray(new_vertices)


def _voronoi_cells(centers, poly):
    """Celdas de Voronoi recortadas a poly. Lista de geometrias.

    Si P es convexo, cada celda es convexa y cubrir sus vertices equivale
    a cubrirla entera.
    """
    import numpy as np
    from scipy.spatial import Voronoi
    from shapely.geometry import Polygon
    C = np.asarray(centers, dtype=float)
    n = len(C)
    # Qhull necesita puntos no degenerados: jitter diminuto si hace falta
    try:
        vor = Voronoi(C)
    except Exception:
        C = C + np.random.RandomState(0).randn(*C.shape) * 1e-9
        vor = Voronoi(C)
    minx, miny, maxx, maxy = poly.bounds
    diag = float(np.hypot(maxx - minx, maxy - miny)) + 1.0
    try:
        regions, verts = _finite_voronoi_regions(vor, diag)
    except Exception:
        return [poly for _ in range(n)]  # fallback conservador
    cells = []
    for reg in regions:
        try:
            cell = Polygon([verts[v] for v in reg])
        except Exception:
            cell = poly
        try:
            inter = cell.intersection(poly)
        except Exception:
            inter = poly
        cells.append(inter)
    return cells


def _cell_constraint_points(cell, max_pts=24):
    """Vertices (+ puntos medios de aristas largas) de la celda."""
    import numpy as np
    if cell.is_empty:
        return np.zeros((0, 2))
    geoms = [cell] if cell.geom_type == "Polygon" else list(cell.geoms)
    pts = []
    for g in geoms:
        if g.is_empty:
            continue
        coords = list(g.exterior.coords)[:-1]
        pts.extend(coords)
        # refuerza aristas largas (no convexo / bordes curvos)
        for (x1, y1), (x2, y2) in zip(coords, coords[1:] + coords[:1]):
            L = float(np.hypot(x2 - x1, y2 - y1))
            if L > 0.1:
                for f in (0.5,):
                    pts.append((x1 + f * (x2 - x1), y1 + f * (y2 - y1)))
    arr = np.array(pts, dtype=float) if pts else np.zeros((0, 2))
    if len(arr) > max_pts:  # diezma uniforme
        idx = np.linspace(0, len(arr) - 1, max_pts).astype(int)
        arr = arr[idx]
    # punto interior por seguridad (celdas no convexas)
    try:
        from shapely.geometry import Point as SPoint
        geoms0 = [cell] if cell.geom_type == "Polygon" \
            else list(cell.geoms)
        for g in geoms0:
            if not g.is_empty:
                p = g.representative_point()
                arr = np.vstack([arr, [[p.x, p.y]]]) if len(arr) else \
                    np.array([[p.x, p.y]])
                break
    except Exception:
        pass
    return arr


def voronoi_refine(poly, R, centers, tour, n_rounds=6, max_step=None,
                   margin=1e-6, verbose=True, verify_tol=1e-6,
                   shape="circle"):
    """Mueve centros en continuo minimizando el tour, cobertura preservada.

    Cada ronda: celdas de Voronoi cap P -> cada centro se mueve (SLSQP)
    para acercarse a su segmento del tour sin soltar ningun vertice de su
    celda -> re-opt 2-opt -> verificacion continua (revierte si rompe).
    shape="square": las cajas son convexas, asi que cubrir los vertices de
    la celda (convexa si P es convexo) equivale a cubrirla entera.
    """
    C = np.asarray(centers, dtype=float)
    t = list(tour)
    if max_step is None:
        max_step = R / 4.0
    L = tour_length(C, t, closed=True)
    hist = [{"round": -1, "k": len(C), "length": L}]
    best = {"centers": C.copy(), "tour": list(t), "length": L, "k": len(C)}
    if verbose:
        print(f"[vor 0] k={len(C)} L={L:.4f}")

    for rd in range(int(n_rounds)):
        cells = _voronoi_cells(C, poly)
        # drops exactos: celda vacia -> centro redundante
        keep = [i for i, c in enumerate(cells) if not c.is_empty]
        if len(keep) < len(C):
            C = C[np.array(keep)]
            t = [keep.index(v) for v in t if v in keep]
            if verbose:
                print(f"[vor {rd}] drops: {len(cells)-len(keep)}")
        arr = C.copy()
        moved = 0
        for pos, node in enumerate(t):
            pts = _cell_constraint_points(cells[node] if node < len(cells)
                                          else poly)
            if len(pts) == 0:
                continue
            prev = arr[t[pos - 1]]
            nxt = arr[t[(pos + 1) % len(t)]]
            new = _best_position(prev, nxt, pts, R - margin, arr[node],
                                 max_step=max_step, shape=shape)
            if float(np.linalg.norm(new - arr[node])) > 1e-9:
                moved += 1
            arr[node] = new
        C = arr
        D = dist_matrix(C)
        t = two_opt(t, D)
        L = tour_length(C, t, closed=True)
        from .verify import check_coverage
        chk = check_coverage(poly, C, R, shape=shape)
        ok = bool(chk["covered"] or chk["uncovered_ratio"] < verify_tol)
        hist.append({"round": rd, "k": len(C), "length": L, "moved": moved,
                     "covered": ok, "uncovered_ratio": chk["uncovered_ratio"]})
        if verbose:
            print(f"[vor {rd+1}] k={len(C)} L={L:.4f} movidos={moved} "
                  f"cubierto={ok} ratio={chk['uncovered_ratio']:.2e}")
        if ok:
            if L < best["length"] - 1e-9:
                best = {"centers": C.copy(), "tour": list(t),
                        "length": L, "k": len(C)}
        else:
            C = best["centers"].copy()
            t = list(best["tour"])
            if verbose:
                print(f"[vor {rd+1}] revierto")
            break
    from .verify import check_coverage as _cc
    chk = _cc(poly, best["centers"], R, shape=shape)
    return {"centers": best["centers"], "tour": best["tour"],
            "length": best["length"], "k": best["k"],
            "history": hist, "check": chk}


# ---------------- backbone (ruta-primero) ----------------

def lawnmower_backbone(poly, spacing):
    """Serpentina horizontal sobre el bbox como LineString."""
    from shapely.geometry import LineString
    minx, miny, maxx, maxy = poly.bounds
    ys = np.arange(miny, maxy + 0.5 * spacing, spacing)
    pts = []
    for k, y in enumerate(ys):
        if k % 2 == 0:
            pts += [(minx, y), (maxx, y)]
        else:
            pts += [(maxx, y), (minx, y)]
    return LineString(pts)


def backbone_cover(poly, R, demands, candidates, lambdas=(0.0, 0.5, 2.0, 8.0),
                   solver="auto", time_limit=120, tsp_time_limit=60,
                   verbose=True, shape="circle"):
    """Barrido conjunto aproximado: coste_j = 1 + lam*dist(j, backbone).

    lam=0 reproduce el dos-fases (min k). lam>0 alinea centros con la
    serpentina -> tours mas cortos a costa de (quiza) mas k.
    Devuelve el mejor por longitud de tour entre los verificados.
    shape: "circle" o "square" (cobertura y verificacion con cuadrados).
    """
    from shapely.geometry import Point
    from .cover import build_coverage
    from .solvers import solve_greedy, solve_set_cover
    from .tsp import solve_tour
    from .verify import check_coverage
    dem = np.asarray(demands, dtype=float)
    cand = np.asarray(candidates, dtype=float)
    cover_lists, infeas = build_coverage(dem, cand, R, shape=shape)
    if infeas:
        return {"feasible": False,
                "message": f"{len(infeas)} demandas sin candidato"}
    bone = lawnmower_backbone(poly, 2 * R)
    d_bone = np.array([bone.distance(Point(p)) for p in cand])
    results = []
    for lam in lambdas:
        w = 1.0 + float(lam) * d_bone
        g = solve_greedy(cover_lists, len(cand), weights=w)
        sol = solve_set_cover(cover_lists, len(cand), method=solver,
                              time_limit=time_limit,
                              warm_start=g["selected"] if g["feasible"]
                              else None, weights=w)
        if not sol["feasible"]:
            continue
        sel = np.array(sol["selected"], dtype=int)
        t = solve_tour(cand[sel], method="auto",
                       time_limit=tsp_time_limit)
        chk = check_coverage(poly, cand[sel], R, shape=shape)
        row = {"lam": lam, "k": len(sel), "length": t["length"],
               "method": t["method"], "covered": chk["covered"],
               "ratio": chk["uncovered_ratio"], "selected": sel,
               "tour_local": t["tour"]}
        results.append(row)
        if verbose:
            print(f"[lam={lam}] k={len(sel)} L={t['length']:.4f} "
                  f"cubierto={chk['covered']} ({t['method']})")
    ok = [r for r in results if r["covered"]]
    if not ok:
        return {"feasible": False, "results": results}
    best = min(ok, key=lambda r: r["length"])
    return {"feasible": True, "results": results, "best": best,
            "centers": cand[np.array(best["selected"], dtype=int)],
            "tour": list(best["tour_local"]), "length": best["length"],
            "k": best["k"]}

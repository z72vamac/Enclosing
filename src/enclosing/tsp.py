"""TSP euclideo sobre los centros: tour cerrado de minima distancia."""
import numpy as np
from scipy.spatial import distance_matrix


def dist_matrix(centers):
    return distance_matrix(np.asarray(centers, dtype=float),
                           np.asarray(centers, dtype=float))


def tour_length(centers, tour, closed=True):
    C = np.asarray(centers, dtype=float)
    t = list(tour)
    if closed:
        t = t + [t[0]]
    L = 0.0
    for a, b in zip(t[:-1], t[1:]):
        L += float(np.linalg.norm(C[a] - C[b]))
    return L


def nearest_neighbor(D, start=0):
    n = D.shape[0]
    unvis = set(range(n))
    unvis.discard(start)
    tour = [start]
    cur = start
    while unvis:
        nxt = min(unvis, key=lambda j: D[cur, j])
        tour.append(nxt)
        unvis.discard(nxt)
        cur = nxt
    return tour


def two_opt(tour, D, max_passes=20):
    """Mejora 2-opt clasica (cerrada)."""
    best = list(tour)
    n = len(best)
    if n < 4:
        return best
    def L(t):
        s = 0.0
        for k in range(n):
            s += D[t[k], t[(k + 1) % n]]
        return s
    best_L = L(best)
    for _ in range(max_passes):
        improved = False
        for i in range(n - 1):
            for j in range(i + 2, n if i > 0 else n - 1):
                a, b = best[i], best[(i + 1) % n]
                c, d = best[j], best[(j + 1) % n]
                gain = (D[a, b] + D[c, d]) - (D[a, c] + D[b, d])
                if gain > 1e-12:
                    best[i + 1:j + 1] = reversed(best[i + 1:j + 1])
                    best_L -= gain
                    improved = True
        if not improved:
            break
    return best


def solve_tsp_heuristic(centers, n_starts=None):
    """Multi-arranque NN + 2-opt. Rapido y buena calidad."""
    C = np.asarray(centers, dtype=float)
    n = len(C)
    D = dist_matrix(C)
    if n_starts is None:
        n_starts = min(n, 10)
    starts = np.linspace(0, n - 1, n_starts, dtype=int)
    best, best_L = None, np.inf
    for s in starts:
        t = two_opt(nearest_neighbor(D, start=int(s)), D)
        Lt = tour_length(C, t, closed=True)
        if Lt < best_L:
            best, best_L = t, Lt
    return {"tour": best, "length": float(best_L), "method": "nn+2opt"}


def solve_tsp_exact_gurobi(centers, time_limit=120, warm_start=None,
                           output_flag=0):
    """TSP exacto MTZ con Gurobi. Solo para n moderado (~<=120)."""
    import gurobipy as gp
    C = np.asarray(centers, dtype=float)
    n = len(C)
    D = dist_matrix(C)
    m = gp.Model()
    m.setParam("OutputFlag", output_flag)
    if time_limit is not None:
        m.setParam("TimeLimit", float(time_limit))
    x = m.addVars(n, n, vtype="B", name="x")
    for i in range(n):
        x[i, i].UB = 0
        m.addConstr(gp.quicksum(x[i, j] for j in range(n)) == 1)
        m.addConstr(gp.quicksum(x[j, i] for j in range(n)) == 1)
    u = m.addVars(n, vtype="C", lb=0, ub=n - 1, name="u")
    m.addConstr(u[0] == 0)
    for i in range(1, n):
        for j in range(1, n):
            if i != j:
                m.addConstr(u[i] - u[j] + n * x[i, j] <= n - 1)
    m.setObjective(gp.quicksum(D[i, j] * x[i, j]
                               for i in range(n) for j in range(n)))
    if warm_start:
        t = list(warm_start) + [warm_start[0]]
        for a, b in zip(t[:-1], t[1:]):
            x[a, b].Start = 1.0
    m.optimize()
    # reconstruye tour desde x
    succ = {}
    try:
        for i in range(n):
            for j in range(n):
                if x[i, j].X > 0.5:
                    succ[i] = j
                    break
        tour = [0]
        for _ in range(n - 1):
            tour.append(succ[tour[-1]])
    except Exception:
        return {"tour": warm_start, "length": None, "method": "gurobi-fail",
                "status": m.status}
    L = tour_length(C, tour, closed=True)
    return {"tour": tour, "length": L, "method": "gurobi-mtz",
            "status": m.status,
            "optimal": bool(m.status == 2 and abs(m.ObjBound - m.ObjVal) < 1e-6)}


def solve_tour(centers, method="auto", time_limit=120):
    """Auto: heuristica y, si n<=120, pulido exacto con warm start."""
    heur = solve_tsp_heuristic(centers)
    if method == "heuristic":
        return heur
    n = len(np.asarray(centers))
    if method in ("auto", "exact") and n <= 120:
        try:
            exact = solve_tsp_exact_gurobi(centers, time_limit=time_limit,
                                           warm_start=heur["tour"])
            if exact.get("length") is not None and exact["length"] <= heur["length"] + 1e-9:
                exact["heuristic_length"] = heur["length"]
                return exact
        except Exception as e:
            pass
    return heur

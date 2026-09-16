"""Solvers para el set cover binario: min sum x s.a. A x >= 1."""
import numpy as np
from scipy.sparse import csr_matrix
from scipy.optimize import milp, LinearConstraint, Bounds


def solve_greedy(cover_lists, Nc, weights=None):
    """Greedy clasico. Con weights: max (nuevos cubiertos)/coste."""
    Nd = len(cover_lists)
    w = np.ones(Nc) if weights is None else np.asarray(weights, dtype=float)
    # invertido: candidatos -> demandas que cubre
    cand_to_dem = [[] for _ in range(Nc)]
    for i, lst in enumerate(cover_lists):
        for j in lst:
            cand_to_dem[j].append(i)
    uncovered = np.ones(Nd, dtype=bool)
    chosen = []
    used = np.zeros(Nc, dtype=bool)
    while uncovered.any():
        best_j, best_gain = -1, -1.0
        for j in range(Nc):
            if used[j]:
                continue
            gain = sum(1 for i in cand_to_dem[j] if uncovered[i]) / max(w[j], 1e-12)
            if gain > best_gain:
                best_gain, best_j = gain, j
        if best_j < 0 or best_gain <= 0:
            break  # infactible
        chosen.append(best_j)
        used[best_j] = True
        for i in cand_to_dem[best_j]:
            uncovered[i] = False
    feasible = not uncovered.any()
    cost = float(sum(w[j] for j in chosen))
    return {"selected": sorted(chosen), "objective": len(chosen),
            "cost": cost, "feasible": feasible, "method": "greedy"}


def _to_csr(cover_lists, Nd, Nc):
    rows, cols, data = [], [], []
    for i, lst in enumerate(cover_lists):
        for j in lst:
            rows.append(i)
            cols.append(j)
            data.append(1.0)
    return csr_matrix((data, (rows, cols)), shape=(Nd, Nc))


def solve_exact_highs(cover_lists, Nc, time_limit=None, warm_start=None,
                      weights=None):
    """Set cover exacto via scipy.optimize.milp (HiGHS incluido)."""
    Nd = len(cover_lists)
    A = _to_csr(cover_lists, Nd, Nc)
    c = np.ones(Nc) if weights is None else np.asarray(weights, dtype=float)
    # A x >= 1  ->  -A x <= -1
    constraints = LinearConstraint(-A, lb=-np.inf, ub=-np.ones(Nd))
    bounds = Bounds(lb=np.zeros(Nc), ub=np.ones(Nc))
    integrality = np.ones(Nc)
    options = {}
    if time_limit is not None:
        options["time_limit"] = float(time_limit)
    res = milp(c=c, constraints=constraints, bounds=bounds,
               integrality=integrality, options=options)
    if res.status not in (0,):
        return {"selected": [], "objective": None, "feasible": False,
                "method": "highs", "status": res.status, "message": str(res.message)}
    x = np.round(res.x).astype(int)
    selected = sorted(int(j) for j in range(Nc) if x[j] > 0)
    return {"selected": selected, "objective": len(selected),
            "feasible": True, "method": "highs", "status": res.status}


def solve_exact_gurobi(cover_lists, Nc, time_limit=None, warm_start=None,
                       output_flag=0, weights=None):
    import gurobipy as gp
    m = gp.Model()
    m.setParam("OutputFlag", output_flag)
    if time_limit is not None:
        m.setParam("TimeLimit", float(time_limit))
    w = np.ones(Nc) if weights is None else np.asarray(weights, dtype=float)
    x = m.addVars(Nc, vtype="B", name="x")
    m.setObjective(gp.quicksum(w[j] * x[j] for j in range(Nc)))
    for i, lst in enumerate(cover_lists):
        if len(lst) == 0:
            return {"selected": [], "objective": None, "feasible": False,
                    "method": "gurobi", "message": f"demanda {i} sin candidatos"}
        m.addConstr(gp.quicksum(x[j] for j in lst) >= 1)
    if warm_start:
        for j in warm_start:
            if 0 <= j < Nc:
                x[j].Start = 1.0
    m.optimize()
    status = m.status
    # 2 = OPTIMAL
    if status != 2:
        return {"selected": [], "objective": None, "feasible": False,
                "method": "gurobi", "status": status}
    selected = sorted(j for j in range(Nc) if x[j].X > 0.5)
    return {"selected": selected, "objective": len(selected),
            "feasible": True, "method": "gurobi", "status": status}


def solve_set_cover(cover_lists, Nc, method="auto", time_limit=120,
                    warm_start=None, weights=None):
    """Despacha solver: gurobi -> highs -> greedy.

    method: 'auto' | 'gurobi' | 'highs' | 'greedy'
    weights: coste por candidato (defecto 1). Con pesos, el goloso usa
    ratio cobertura/coste y los exactos minimizan coste total.
    """
    if method in ("auto", "gurobi"):
        try:
            r = solve_exact_gurobi(cover_lists, Nc, time_limit=time_limit,
                                   warm_start=warm_start, weights=weights)
            if r["feasible"]:
                return r
            if method == "gurobi":
                return r
        except Exception as e:  # sin licencia / sin gurobipy
            if method == "gurobi":
                return {"selected": [], "objective": None, "feasible": False,
                        "method": "gurobi", "message": str(e)}
    if method in ("auto", "highs"):
        try:
            r = solve_exact_highs(cover_lists, Nc, time_limit=time_limit,
                                  weights=weights)
            if r["feasible"]:
                return r
            if method == "highs":
                return r
        except Exception as e:
            if method == "highs":
                return {"selected": [], "objective": None, "feasible": False,
                        "method": "highs", "message": str(e)}
    g = solve_greedy(cover_lists, Nc, weights=weights)
    return g

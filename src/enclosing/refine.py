"""Bucle iterativo: resolver set cover discreto -> verificar continuo -> anadir demanda."""
import numpy as np
from .sampling import demand_points, candidate_centers, sample_points_in_geom
from .cover import build_coverage
from .solvers import solve_greedy, solve_set_cover
from .verify import check_coverage
from .geometry import lower_bounds


def _dedupe(points, tol):
    if len(points) == 0:
        return points
    key = np.round(np.asarray(points) / (tol + 1e-12)).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    return np.asarray(points)[np.sort(idx)]


def iterative_cover(poly, R, h0=None, s0=None, max_iter=8,
                    solver="auto", time_limit=120, verbose=True,
                    buffer_resolution=64, cand_method="hex"):
    """Cubre poly con discos de radio R.

    h0: paso de demanda inicial (defecto R/2). s0: paso de candidatos (defecto R/2).
    Devuelve dict con centers, k, demands, candidates, bounds, history, check.
    """
    if h0 is None:
        h0 = R / 2.0
    if s0 is None:
        s0 = R / 2.0
    bounds = lower_bounds(poly, R)
    s = float(s0)
    h_new = float(h0)
    candidates = candidate_centers(poly, R, s, method=cand_method)
    demands = demand_points(poly, h0)
    history = []
    result = None

    for it in range(int(max_iter)):
        cover_lists, uncovered_idx = build_coverage(demands, candidates, R)
        if uncovered_idx:
            # candidatos demasiado gruesos: densificar y reintentar
            if s < R / 8.0:
                history.append({"iter": it, "status": "infeasible_discrete",
                                "n_dem": len(demands), "n_cand": len(candidates)})
                break
            s = s / 2.0
            candidates = candidate_centers(poly, R, s, method=cand_method)
            if verbose:
                print(f"[iter {it}] {len(uncovered_idx)} demandas sin candidato: "
                      f"densifico candidatos a s={s:.4f} ({len(candidates)}).")
            continue
        greedy = solve_greedy(cover_lists, len(candidates))
        sol = solve_set_cover(cover_lists, len(candidates), method=solver,
                              time_limit=time_limit,
                              warm_start=greedy["selected"] if greedy["feasible"] else None)
        if not sol["feasible"]:
            history.append({"iter": it, "status": "solver_infeasible",
                            "n_dem": len(demands), "n_cand": len(candidates)})
            break
        centers = candidates[np.array(sol["selected"], dtype=int)]
        chk = check_coverage(poly, centers, R, buffer_resolution=buffer_resolution)
        history.append({"iter": it, "k": len(centers), "n_dem": len(demands),
                        "n_cand": len(candidates), "method": sol["method"],
                        "uncovered_area": chk["uncovered_area"],
                        "uncovered_ratio": chk["uncovered_ratio"],
                        "covered": chk["covered"]})
        if verbose:
            print(f"[iter {it}] k={len(centers)} dem={len(demands)} "
                  f"cand={len(candidates)} uncovered={chk['uncovered_ratio']:.3e} "
                  f"({sol['method']})")
        if chk["covered"]:
            result = {"centers": centers, "k": len(centers), "demands": demands,
                      "candidates": candidates, "selected": sol["selected"],
                      "bounds": bounds, "history": history, "check": chk,
                      "cover_lists": cover_lists}
            break
        # anadir demanda donde falta cobertura
        h_new = h_new / 2.0
        new_pts = sample_points_in_geom(chk["uncovered_geom"], h_new)
        # filtra los ya presentes
        demands = _dedupe(np.vstack([demands, new_pts]) if len(new_pts) else demands,
                          tol=h_new / 4.0)
        result = {"centers": centers, "k": len(centers), "demands": demands,
                  "candidates": candidates, "selected": sol["selected"],
                  "bounds": bounds, "history": history, "check": chk,
                  "cover_lists": cover_lists}
    return result

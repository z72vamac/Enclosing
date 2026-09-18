"""Densidad de demanda/candidatos en el irregular R=1.0: A base, B dem densa, C todo denso."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import time
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.sampling import demand_points, candidate_centers
from enclosing.cover import build_coverage
from enclosing.solvers import solve_set_cover, solve_greedy
from enclosing.tsp import solve_tour
from enclosing.joint import solve_joint_dfj
from enclosing.verify import check_coverage

poly = make_polygon([(0, 0), (5, 0), (6, 3), (3, 5), (1, 4.5), (0, 2.5)])
R = 1.0

for name, h, s in (("A base (h=0.5,s=1.0)", 0.5, 1.0),
                   ("B dem densa (h=0.25,s=1.0)", 0.25, 1.0),
                   ("C todo denso (h=0.25,s=0.5)", 0.25, 0.5)):
    dem = demand_points(poly, h)
    cand = candidate_centers(poly, R, s)
    cl, ui = build_coverage(dem, cand, R)
    assert not ui, f"{name}: {len(ui)} demandas sin candidato"
    depot = int(np.argmin(np.linalg.norm(cand - np.array([2.8, 2.2]), axis=1)))
    g = solve_greedy(cl, len(cand))
    cov = solve_set_cover(cl, len(cand), method="auto", time_limit=120,
                          warm_start=g["selected"])
    sel2 = sorted(set(cov["selected"]) | {depot})
    t2 = solve_tour(cand[np.array(sel2)], method="auto", time_limit=120)
    warm_global = [sel2[i] for i in t2["tour"]]
    t0 = time.time()
    j = solve_joint_dfj(cand, cl, depot_idx=depot, time_limit=300,
                        warm_y=sel2, warm_tour=warm_global)
    dt = time.time() - t0
    chk = check_coverage(poly, cand[np.array(j["selected"])], R,
                         buffer_resolution=256) if j["feasible"] else None
    print(f"{name}: dem={len(dem)} cand={len(cand)} "
          f"2fases k={len(sel2)} L={t2['length']:.4f} | "
          f"joint k={j.get('k')} L={round(j.get('length', float('nan')), 4)} "
          f"gap={j.get('gap')} optimo={j.get('optimal')} t={dt:.0f}s | "
          f"ratio256={chk['uncovered_ratio']:.2e}" if chk else
          f"{name}: dem={len(dem)} cand={len(cand)} NO FACTIBLE")

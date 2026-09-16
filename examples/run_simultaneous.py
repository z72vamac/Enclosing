"""Simultaneo DFJ vs dos-fases en instancia media (cuadrado 4x4, R=1)."""
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

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 1.0
dem = demand_points(poly, 0.5)
cand = candidate_centers(poly, R, 1.0)
print("dem:", len(dem), "cand:", len(cand))
cl, ui = build_coverage(dem, cand, R)
assert not ui

centroid = np.array([2.0, 2.0])
depot = int(np.argmin(np.linalg.norm(cand - centroid, axis=1)))

# dos-fases con depot
g = solve_greedy(cl, len(cand))
cov = solve_set_cover(cl, len(cand), method="auto", time_limit=60,
                      warm_start=g["selected"])
sel2 = sorted(set(cov["selected"]) | {depot})
t2 = solve_tour(cand[np.array(sel2)], method="auto", time_limit=120)
print(f"dos-fases: k={len(sel2)} L={t2['length']:.4f}")

# simultaneo DFJ (warm start = dos-fases)
th = t2["tour"]
warm_global = [sel2[i] for i in th]
t0 = time.time()
j = solve_joint_dfj(cand, cl, depot_idx=depot, time_limit=300,
                    warm_y=sel2, warm_tour=warm_global)
dt = time.time() - t0
print(f"joint DFJ: k={j.get('k')} L={round(j.get('length', float('nan')), 4)} "
      f"gap={j.get('gap')} optimo={j.get('optimal')} t={dt:.1f}s")
if j["feasible"]:
    chk = check_coverage(poly, cand[np.array(j["selected"])], R)
    print("cubierto:", chk["covered"], f"{chk['uncovered_ratio']:.2e}")
    print("mejora tour:", round(t2["length"] - j["length"], 4))

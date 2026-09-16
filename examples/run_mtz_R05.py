"""MTZ en R=0.5 con mismo warm-start y limite que el DFJ."""
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
from enclosing.joint import solve_joint_mtz

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 0.5
dem = demand_points(poly, R / 2)
cand = candidate_centers(poly, R, R)
cl, ui = build_coverage(dem, cand, R)
assert not ui
depot = int(np.argmin(np.linalg.norm(cand - np.array([2.0, 2.0]), axis=1)))

g = solve_greedy(cl, len(cand))
cov = solve_set_cover(cl, len(cand), method="auto", time_limit=120,
                      warm_start=g["selected"])
sel2 = sorted(set(cov["selected"]) | {depot})
t2 = solve_tour(cand[np.array(sel2)], method="auto", time_limit=120)
print(f"dos-fases: k={len(sel2)} L={t2['length']:.4f}")

warm_global = [sel2[i] for i in t2["tour"]]
t0 = time.time()
j = solve_joint_mtz(cand, cl, depot_idx=depot, time_limit=600,
                    warm_y=sel2, warm_tour=warm_global)
print(f"MTZ: factible={j['feasible']} k={j.get('k')} "
      f"L={round(j.get('length', float('nan')), 4) if j['feasible'] else None} "
      f"gap={j.get('gap')} optimo={j.get('optimal')} t={time.time()-t0:.1f}s")

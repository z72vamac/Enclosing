"""Joint exacto vs dos-fases en instancia pequena (cuadrado 2x2, R=1)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.sampling import demand_points, candidate_centers
from enclosing.cover import build_coverage
from enclosing.solvers import solve_set_cover, solve_greedy
from enclosing.tsp import solve_tour
from enclosing.joint import solve_joint_exact

poly = make_polygon([(0, 0), (2, 0), (2, 2), (0, 2)])
R = 1.0
dem = demand_points(poly, 0.5)
cand = candidate_centers(poly, R, 1.0)
print("dem:", len(dem), "cand:", len(cand))
cl, ui = build_coverage(dem, cand, R)
assert not ui

# base = centroide -> depot del dron
centroid = np.array([1.0, 1.0])
depot = int(np.argmin(np.linalg.norm(cand - centroid, axis=1)))
print("depot:", depot, cand[depot].tolist())

# --- dos fases (tour obligado a pasar por depot) ---
g = solve_greedy(cl, len(cand))
cov = solve_set_cover(cl, len(cand), method="auto", time_limit=60,
                      warm_start=g["selected"])
sel2 = sorted(set(cov["selected"]) | {depot})
t2 = solve_tour(cand[np.array(sel2)], method="heuristic")
# remapea a indices globales
inv = {v: k for k, v in enumerate(sel2)}
L2 = t2["length"]
print(f"dos-fases: k={len(sel2)} L={L2:.4f} ({cov['method']})")

# --- conjunto exacto ---
warm_tour_local = [inv[v] for v in [depot]]  # no usado; warm via heur global
# warm start: tour heuristico remapeado a global
th = solve_tour(cand[np.array(sel2)], method="heuristic")["tour"]
warm_global = [sel2[i] for i in th]
j = solve_joint_exact(cand, cl, depot_idx=depot, time_limit=120,
                      warm_y=sel2, warm_tour=warm_global)
print("joint:", j["k"], round(j["length"], 4), "optimo:", j.get("optimal"),
      "factible:", j["feasible"])
print("mejora:", round(L2 - j["length"], 4) if j["feasible"] else None)

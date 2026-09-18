"""Irregular 2 fases vs conjunto (hexagono irregular, R=1.0, circulos)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from enclosing.geometry import make_polygon
from enclosing.sampling import demand_points, candidate_centers
from enclosing.cover import build_coverage
from enclosing.solvers import solve_set_cover, solve_greedy
from enclosing.tsp import solve_tour
from enclosing.joint import solve_joint_dfj
from enclosing.verify import check_coverage

poly = make_polygon([(0, 0), (5, 0), (6, 3), (3, 5), (1, 4.5), (0, 2.5)])
R = 1.0
dem = demand_points(poly, 0.5)
cand = candidate_centers(poly, R, 1.0)
print("dem:", len(dem), "cand:", len(cand))
cl, ui = build_coverage(dem, cand, R)
assert not ui, f"{len(ui)} demandas sin candidato"
centroid = np.array([2.8, 2.2])
depot = int(np.argmin(np.linalg.norm(cand - centroid, axis=1)))
print("depot:", depot, cand[depot].tolist())

# --- dos fases ---
g = solve_greedy(cl, len(cand))
cov = solve_set_cover(cl, len(cand), method="auto", time_limit=120,
                      warm_start=g["selected"])
sel2 = sorted(set(cov["selected"]) | {depot})
t2 = solve_tour(cand[np.array(sel2)], method="auto", time_limit=120)
tour2 = [sel2[i] for i in t2["tour"]]
C2 = cand[np.array(sel2)]
order2 = [sel2.index(v) for v in tour2]
chk2 = check_coverage(poly, C2, R)
print(f"dos-fases: k={len(sel2)} L={t2['length']:.4f} ({t2['method']}) "
      f"cubierto={chk2['covered']} ratio={chk2['uncovered_ratio']:.2e}")

# --- conjunto DFJ ---
warm_global = [sel2[i] for i in t2["tour"]]
t0 = time.time()
j = solve_joint_dfj(cand, cl, depot_idx=depot, time_limit=300,
                    warm_y=sel2, warm_tour=warm_global)
dt = time.time() - t0
print(f"joint DFJ: factible={j['feasible']} k={j.get('k')} "
      f"L={round(j.get('length', float('nan')), 4)} gap={j.get('gap')} "
      f"optimo={j.get('optimal')} t={dt:.1f}s")
Cj = cand[np.array(j["selected"])]
orderj = [j["selected"].index(v) for v in j["tour"]]
chkJ = check_coverage(poly, Cj, R)
print(f"joint: cubierto={chkJ['covered']} ratio={chkJ['uncovered_ratio']:.2e}")
print(f"mejora tour dos-fases -> conjunto: {t2['length'] - j['length']:.4f} "
      f"(dk={j['k'] - len(sel2):+d})")

# --- figura ---
fig, ax = plt.subplots(1, 2, figsize=(12, 5.5))
for a, C, order, ttl in ((ax[0], C2, order2,
                          f"Dos fases k={len(C2)} L={t2['length']:.2f}"),
                         (ax[1], Cj, orderj,
                          f"Conjunto k={len(Cj)} L={j['length']:.2f}")):
    a.set_aspect("equal")
    x, y = poly.exterior.xy
    a.fill(x, y, alpha=0.12)
    a.plot(x, y, "k-", lw=1.2)
    for (cx, cy) in C:
        a.add_patch(Circle((cx, cy), R, alpha=0.18, edgecolor="C0",
                           facecolor="C0"))
    tt = list(order) + [order[0]]
    a.plot(C[tt, 0], C[tt, 1], "C1-", lw=1.0, label="tour")
    a.plot(C[:, 0], C[:, 1], "C0o", ms=3)
    a.set_title(ttl)
    a.legend(loc="best")
plt.tight_layout()
out = str(Path(__file__).resolve().parent / "irregular_joint_R1.png")
fig.savefig(out, dpi=150)
print("figura:", out)

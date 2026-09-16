"""Figura comparativa dos-fases vs simultaneo (cuadrado 4x4, R=1)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from enclosing.geometry import make_polygon
from enclosing.sampling import demand_points, candidate_centers
from enclosing.cover import build_coverage
from enclosing.solvers import solve_set_cover, solve_greedy
from enclosing.tsp import solve_tour
from enclosing.joint import solve_joint_dfj

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 1.0
dem = demand_points(poly, 0.5)
cand = candidate_centers(poly, R, 1.0)
cl, ui = build_coverage(dem, cand, R)
assert not ui
depot = int(np.argmin(np.linalg.norm(cand - np.array([2.0, 2.0]), axis=1)))

g = solve_greedy(cl, len(cand))
cov = solve_set_cover(cl, len(cand), method="auto", time_limit=60,
                      warm_start=g["selected"])
sel2 = sorted(set(cov["selected"]) | {depot})
t2 = solve_tour(cand[np.array(sel2)], method="auto", time_limit=120)
tour2 = [sel2[i] for i in t2["tour"]]              # a indices globales
C2 = cand[np.array(sel2)]
order2 = [sel2.index(v) for v in tour2]           # de vuelta a locales

warm_global = [sel2[i] for i in t2["tour"]]
j = solve_joint_dfj(cand, cl, depot_idx=depot, time_limit=300,
                    warm_y=sel2, warm_tour=warm_global)
Cj = cand[np.array(j["selected"])]
orderj = [j["selected"].index(v) for v in j["tour"]]  # idem

fig, ax = plt.subplots(1, 2, figsize=(11, 5.5))
for a, C, order, ttl in ((ax[0], C2, order2,
                          f"Dos fases k={len(C2)} L={t2['length']:.2f}"),
                         (ax[1], Cj, orderj,
                          f"Simultaneo k={len(Cj)} L={j['length']:.2f}")):
    a.set_aspect("equal")
    x, y = poly.exterior.xy
    a.fill(x, y, alpha=0.12)
    a.plot(x, y, "k-", lw=1.2)
    for (cx, cy) in C:
        a.add_patch(Circle((cx, cy), R, alpha=0.18, edgecolor="C0",
                           facecolor="C0"))
    tt = list(order) + [order[0]]
    a.plot(C[tt, 0], C[tt, 1], "C1-", lw=1.2, label="tour")
    a.plot(C[:, 0], C[:, 1], "C0o", ms=4)
    a.set_title(ttl)
    a.legend(loc="best")
plt.tight_layout()
out = str(Path(__file__).resolve().parent / "simultaneous_4x4.png")
fig.savefig(out, dpi=150)
print("figura:", out)

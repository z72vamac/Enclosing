"""Pulido del caso denso C (h=0.25,s=0.5) hasta certificado continuo."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.sampling import demand_points, candidate_centers, sample_points_in_geom
from enclosing.cover import build_coverage
from enclosing.solvers import solve_set_cover, solve_greedy
from enclosing.tsp import solve_tour, tour_length
from enclosing.joint import solve_joint_dfj
from enclosing.verify import check_coverage
from enclosing.visualize import plot_tour

poly = make_polygon([(0, 0), (5, 0), (6, 3), (3, 5), (1, 4.5), (0, 2.5)])
R = 1.0
DEM = demand_points(poly, 0.25)
CAND = candidate_centers(poly, R, 0.5)
CL, ui = build_coverage(DEM, CAND, R)
assert not ui
depot = int(np.argmin(np.linalg.norm(CAND - np.array([2.8, 2.2]), axis=1)))
g = solve_greedy(CL, len(CAND))
cov = solve_set_cover(CL, len(CAND), method="auto", time_limit=120,
                      warm_start=g["selected"])
sel_prev = sorted(set(cov["selected"]) | {depot})
t2 = solve_tour(CAND[np.array(sel_prev)], method="auto", time_limit=120)
tour_prev = [sel_prev[i] for i in t2["tour"]]
print(f"2fases denso: k={len(sel_prev)} L={t2['length']:.4f}")

for rd in range(4):
    CL, ui = build_coverage(DEM, CAND, R)
    print(f"[denso {rd}] demandas={len(DEM)} sin-candidato={len(ui)}")
    assert not ui
    if set(tour_prev) != set(sel_prev) or len(tour_prev) != len(sel_prev):
        print(f"[denso {rd}] warm-tour parcial: reconstruyo orden heuristico")
        th = solve_tour(CAND[np.array(sel_prev)], method="heuristic")["tour"]
        tour_prev = [sel_prev[i] for i in th]
    jj = solve_joint_dfj(CAND, CL, depot_idx=depot, time_limit=300,
                         warm_y=sel_prev, warm_tour=tour_prev)
    if not jj.get("feasible"):
        print(f"[denso {rd}] joint infactible: {jj}")
        break
    sel_prev = sorted(jj["selected"])
    tour_prev = list(jj["tour"])
    chk = check_coverage(poly, CAND[np.array(sel_prev)], R,
                         buffer_resolution=256)
    print(f"[denso {rd}] joint: k={jj['k']} L={jj['length']:.4f} "
          f"gap={jj.get('gap')} optimo={jj.get('optimal')} "
          f"ratio256={chk['uncovered_ratio']:.2e} cubierto={chk['covered']}")
    if chk["covered"] or chk["uncovered_ratio"] < 1e-7:
        break
    DEM = np.vstack([DEM, sample_points_in_geom(chk["uncovered_geom"], R / 4)])

Cfin = CAND[np.array(sel_prev)]
inv = {gg: i for i, gg in enumerate(sel_prev)}
tour_loc = [inv[v] for v in tour_prev if v in inv]
Lbest = tour_length(Cfin, tour_loc, closed=True) if len(tour_loc) == len(sel_prev) else float("nan")
plot_tour(poly, Cfin, tour_loc, R=R,
          save=str(Path(__file__).resolve().parent / "irregular_dense_polished.png"),
          title=f"Irregular denso pulido k={len(sel_prev)} L={Lbest:.2f}")
print("figura: irregular_dense_polished.png")

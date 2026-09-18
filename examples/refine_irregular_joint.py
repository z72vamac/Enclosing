"""Refino continuo (Voronoi) de la solucion conjunta en el irregular."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.sampling import demand_points, candidate_centers
from enclosing.cover import build_coverage
from enclosing.solvers import solve_set_cover, solve_greedy
from enclosing.tsp import solve_tour
from enclosing.joint import solve_joint_dfj, voronoi_refine, alternating_refine
from enclosing.verify import check_coverage
from enclosing.visualize import plot_tour

poly = make_polygon([(0, 0), (5, 0), (6, 3), (3, 5), (1, 4.5), (0, 2.5)])
R = 1.0
dem = demand_points(poly, 0.5)
cand = candidate_centers(poly, R, 1.0)
cl, ui = build_coverage(dem, cand, R)
assert not ui
depot = int(np.argmin(np.linalg.norm(cand - np.array([2.8, 2.2]), axis=1)))

g = solve_greedy(cl, len(cand))
cov = solve_set_cover(cl, len(cand), method="auto", time_limit=120,
                      warm_start=g["selected"])
sel2 = sorted(set(cov["selected"]) | {depot})
t2 = solve_tour(cand[np.array(sel2)], method="auto", time_limit=120)
warm_global = [sel2[i] for i in t2["tour"]]
j = solve_joint_dfj(cand, cl, depot_idx=depot, time_limit=300,
                    warm_y=sel2, warm_tour=warm_global)
print(f"joint DFJ: k={j['k']} L={j['length']:.4f} optimo={j.get('optimal')}")

C0 = cand[np.array(j["selected"])]
inv = {g: i for i, g in enumerate(j["selected"])}
tour_local = [inv[v] for v in j["tour"]]
out = voronoi_refine(poly, R, C0, tour_local, n_rounds=6,
                     max_step=R / 4, verbose=True)
print(f"voronoi: k={out['k']} L={out['length']:.4f} "
      f"cubierto={out['check']['covered']} "
      f"ratio={out['check']['uncovered_ratio']:.2e}")
print(f"mejora: {j['length'] - out['length']:.4f}")
plot_tour(poly, out["centers"], out["tour"], R=R,
          save=str(Path(__file__).resolve().parent / "irregular_joint_refined.png"),
          title=f"Irregular refinado k={out['k']} L={out['length']:.2f}")
print("figura: irregular_joint_refined.png")

# --- alternante con guarda densa + chequeo en alta resolucion ---
for reso in (64, 128, 256):
    chk = check_coverage(poly, C0, R, buffer_resolution=reso)
    print(f"joint alta-res {reso}: cubierto={chk['covered']} "
          f"ratio={chk['uncovered_ratio']:.2e}")
guard = demand_points(poly, R / 4)
print("guarda densa:", len(guard))
alt = alternating_refine(poly, R, C0, tour_local, guard, n_rounds=10,
                         try_drops=False, verbose=True, max_step=R / 8)
for reso in (64, 256):
    chk = check_coverage(poly, alt["centers"], R, buffer_resolution=reso)
    print(f"alternante alta-res {reso}: cubierto={chk['covered']} "
          f"ratio={chk['uncovered_ratio']:.2e}")
print(f"alternante: k={alt['k']} L={alt['length']:.4f} "
      f"mejora vs DFJ: {j['length'] - alt['length']:.4f}")
plot_tour(poly, alt["centers"], alt["tour"], R=R,
          save=str(Path(__file__).resolve().parent / "irregular_joint_alt.png"),
          title=f"Irregular alternante k={alt['k']} L={alt['length']:.2f}")
print("figura: irregular_joint_alt.png")

# --- (c) enriquece demanda en las rendijas y re-resuelve el conjunto ---
from enclosing.sampling import sample_points_in_geom
DEM = dem
CAND = cand
CL = cl
sel_prev, tour_prev, L_prev = sel2, warm_global, t2["length"]
for rd in range(3):
    Csel = CAND[np.array(sel_prev)]
    chk = check_coverage(poly, Csel, R, buffer_resolution=256)
    print(f"[pulido {rd}] k={len(sel_prev)} ratio={chk['uncovered_ratio']:.2e}")
    if chk["covered"] or chk["uncovered_ratio"] < 1e-7:
        break
    new_pts = sample_points_in_geom(chk["uncovered_geom"], R / 4)
    DEM = np.vstack([DEM, new_pts])
    CL, ui = build_coverage(DEM, CAND, R)
    assert not ui, "candidatos insuficientes para las rendijas"
    jj = solve_joint_dfj(CAND, CL, depot_idx=depot, time_limit=300,
                         warm_y=sel_prev, warm_tour=tour_prev)
    sel_prev = sorted(jj["selected"])
    tour_prev = list(jj["tour"])
    L_prev = jj["length"]
    print(f"[pulido {rd}] joint: k={jj['k']} L={jj['length']:.4f} "
          f"gap={jj.get('gap')} optimo={jj.get('optimal')}")
Cfin = CAND[np.array(sel_prev)]
inv = {g: i for i, g in enumerate(sel_prev)}
tour_fin = [inv[v] for v in tour_prev]
chk_fin = check_coverage(poly, Cfin, R, buffer_resolution=256)
print(f"FINAL pulido: k={len(sel_prev)} L={L_prev:.4f} "
      f"cubierto={chk_fin['covered']} ratio={chk_fin['uncovered_ratio']:.2e}")
plot_tour(poly, Cfin, tour_fin, R=R,
          save=str(Path(__file__).resolve().parent / "irregular_joint_polished.png"),
          title=f"Irregular pulido k={len(sel_prev)} L={L_prev:.2f}")
print("figura: irregular_joint_polished.png")

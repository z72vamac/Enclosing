"""Alternante sobre la solucion R=0.25 (110 centros)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.refine import iterative_cover
from enclosing.sampling import demand_points
from enclosing.tsp import solve_tour
from enclosing.joint import alternating_refine
from enclosing.visualize import plot_tour

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 0.25
res = iterative_cover(poly, R, h0=R / 2, s0=R, max_iter=1, verbose=False,
                      time_limit=120, buffer_resolution=64)
C0 = np.asarray(res["centers"])
t0 = solve_tour(C0, method="auto", time_limit=120)
print("dos-fases: k =", len(C0), "L =", round(t0["length"], 4))

guard = demand_points(poly, R / 4)  # malla densa: los movimientos la respetan
print("guarda:", len(guard))
out = alternating_refine(poly, R, C0, t0["tour"], guard,
                         n_rounds=15, try_drops=False, verbose=True,
                         max_step=R / 32)
print("alternante: k =", out["k"], "L =", round(out["length"], 4),
      "cubierto:", out["check"]["covered"],
      "ratio:", f"{out['check']['uncovered_ratio']:.2e}")
np.save(str(Path(__file__).resolve().parent / "tour_joint_R025.npy"),
        np.array(out["tour"], dtype=int))
np.save(str(Path(__file__).resolve().parent / "centers_joint_R025.npy"),
        out["centers"])
plot_tour(poly, out["centers"], out["tour"], R=R,
          save=str(Path(__file__).resolve().parent / "tour_joint_R025.png"),
          title=f"Joint A k={out['k']} L={out['length']:.2f}")

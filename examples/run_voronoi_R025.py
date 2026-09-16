"""Voronoi continuo sobre la solucion R=0.25."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.refine import iterative_cover
from enclosing.tsp import solve_tour
from enclosing.joint import voronoi_refine
from enclosing.visualize import plot_tour

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 0.25
res = iterative_cover(poly, R, h0=R / 2, s0=R, max_iter=1, verbose=False,
                      time_limit=120, buffer_resolution=64)
C0 = np.asarray(res["centers"])
t0 = solve_tour(C0, method="auto", time_limit=120)
print("dos-fases: k =", len(C0), "L =", round(t0["length"], 4))

out = voronoi_refine(poly, R, C0, t0["tour"], n_rounds=6,
                     max_step=R / 4, verbose=True)
print("voronoi: k =", out["k"], "L =", round(out["length"], 4),
      "cubierto:", out["check"]["covered"],
      "ratio:", f"{out['check']['uncovered_ratio']:.2e}")
np.save(str(Path(__file__).resolve().parent / "centers_voronoi_R025.npy"),
        out["centers"])
np.save(str(Path(__file__).resolve().parent / "tour_voronoi_R025.npy"),
        np.array(out["tour"], dtype=int))
plot_tour(poly, out["centers"], out["tour"], R=R,
          save=str(Path(__file__).resolve().parent / "tour_voronoi_R025.png"),
          title=f"Voronoi k={out['k']} L={out['length']:.2f}")
print("figura: tour_voronoi_R025.png")

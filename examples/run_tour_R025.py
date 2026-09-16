"""Tour de minima distancia sobre los centros del recubrimiento R=0.25."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.refine import iterative_cover
from enclosing.tsp import solve_tour, tour_length
from enclosing.visualize import plot_tour

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 0.25
res = iterative_cover(poly, R, h0=R / 2, s0=R, max_iter=1, verbose=True,
                      time_limit=120, buffer_resolution=64)
C = np.asarray(res["centers"])
print("k =", len(C))

heur = solve_tour(C, method="heuristic")
print("heuristica:", round(heur["length"], 4))

exact = solve_tour(C, method="auto", time_limit=120)
print("final:", exact["method"], round(exact["length"], 4))
np.save(str(Path(__file__).resolve().parent / "tour_R025.npy"),
        np.array(exact["tour"], dtype=int))

out = str(Path(__file__).resolve().parent / "tour_R025.png")
plot_tour(poly, C, exact["tour"], R=R, save=out,
          title=f"Tour {len(C)} centros R=0.25 L={exact['length']:.2f}")
print("figura:", out)

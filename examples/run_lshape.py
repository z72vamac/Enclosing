"""Ejemplo: forma en L con R=0.8."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from enclosing.geometry import make_polygon
from enclosing.refine import iterative_cover
from enclosing.visualize import plot_solution

poly = make_polygon([(0, 0), (3, 0), (3, 1), (1, 1), (1, 3), (0, 3)])
R = 0.8
res = iterative_cover(poly, R, h0=0.4, s0=0.4, max_iter=8, verbose=True)
print("cotas:", res["bounds"])
print("k =", res["k"], "cubierto:", res["check"]["covered"])
out = Path(__file__).resolve().parent / "lshape.png"
plot_solution(poly, res["centers"], R, demands=res["demands"], save=str(out),
              title=f"L R=0.8 k={res['k']}")
print("figura:", out)

"""Ejemplo: cuadrado 4x4 con R=1."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from enclosing.geometry import make_polygon
from enclosing.refine import iterative_cover
from enclosing.visualize import plot_solution

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 1.0
res = iterative_cover(poly, R, h0=0.5, s0=0.5, max_iter=8, verbose=True)
print("cotas:", res["bounds"])
print("k =", res["k"], "cubierto:", res["check"]["covered"])
out = Path(__file__).resolve().parent / "square.png"
plot_solution(poly, res["centers"], R, demands=res["demands"], save=str(out),
              title=f"Cuadrado 4x4 R=1 k={res['k']}")
print("figura:", out)

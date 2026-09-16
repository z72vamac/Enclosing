"""Cuadrado 4x4 con R=0.25."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from enclosing.geometry import make_polygon
from enclosing.refine import iterative_cover
from enclosing.visualize import plot_solution

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 0.25
res = iterative_cover(poly, R, h0=R / 2, s0=R, max_iter=2, verbose=True,
                      time_limit=120, buffer_resolution=64)
print("cotas:", res["bounds"])
print("k =", res["k"], "cubierto:", res["check"]["covered"])
print("historia:", res["history"])
out = Path(__file__).resolve().parent / "square_R025.png"
plot_solution(poly, res["centers"], R, demands=None, save=str(out),
              title=f"Cuadrado 4x4 R=0.25 k={res['k']}")
print("figura:", out)
print("primeros centros:", res["centers"][:10].tolist())

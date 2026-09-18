"""Cuadrado 4x4 con cuadrados de semilado s=0.75 (orientacion fija)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from enclosing.geometry import make_polygon
from enclosing.refine import iterative_cover
from enclosing.visualize import plot_solution

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 0.75
res = iterative_cover(poly, R, h0=R / 2, s0=R, max_iter=4, verbose=True,
                      time_limit=120, shape="square")
print("cotas:", res["bounds"])
print("k =", res["k"], "cubierto:", res["check"]["covered"],
      "ratio:", res["check"]["uncovered_ratio"])
print("historia:", [(h.get("iter"), h.get("k"), round(h.get("uncovered_ratio", -1), 6))
                    for h in res["history"]])
out = Path(__file__).resolve().parent / "square_S075.png"
plot_solution(poly, res["centers"], R, demands=None, save=str(out),
              title=f"Cuadrado 4x4 s=0.75 k={res['k']}", shape="square")
print("figura:", out)
print("centros:", res["centers"].tolist())

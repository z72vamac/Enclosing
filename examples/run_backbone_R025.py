"""Backbone conjunto en R=0.25: barrido en lam."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.sampling import demand_points, candidate_centers
from enclosing.joint import backbone_cover
from enclosing.visualize import plot_tour

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 0.25
dem = demand_points(poly, R / 2)
cand = candidate_centers(poly, R, R)
print("dem:", len(dem), "cand:", len(cand))
out = backbone_cover(poly, R, dem, cand,
                     lambdas=(0.0, 2.0, 8.0, 30.0),
                     solver="auto", time_limit=120, tsp_time_limit=60,
                     verbose=True)
print("mejor:", out["best"])
np.save(str(Path(__file__).resolve().parent / "centers_backbone_R025.npy"),
        out["centers"])
np.save(str(Path(__file__).resolve().parent / "tour_backbone_R025.npy"),
        np.array(out["tour"], dtype=int))
plot_tour(poly, out["centers"], out["tour"], R=R,
          save=str(Path(__file__).resolve().parent / "tour_backbone_R025.png"),
          title=f"Backbone k={out['k']} L={out['length']:.2f}")

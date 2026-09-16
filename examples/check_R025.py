import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
from enclosing.geometry import make_polygon
from enclosing.refine import iterative_cover
from enclosing.verify import check_coverage

poly = make_polygon([(0, 0), (4, 0), (4, 4), (0, 4)])
R = 0.25
res = iterative_cover(poly, R, h0=R / 2, s0=R, max_iter=1, verbose=False,
                      time_limit=120, buffer_resolution=32)
for reso in (32, 64, 128, 256):
    chk = check_coverage(poly, res["centers"], R, buffer_resolution=reso)
    print(reso, chk["covered"], f"{chk['uncovered_ratio']:.3e}")

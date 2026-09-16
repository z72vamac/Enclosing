from triangle_enclosing import triangle_enclosing
from triangle_enclosing_discrete import triangle_enclosing_discrete

import numpy as np

np.random.seed(2)

A = np.random.uniform(0, 10, (10, 2))
C = A

triangle_enclosing(A, C, L = 40)
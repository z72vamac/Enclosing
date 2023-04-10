from triangle_enclosing import triangle_enclosing

import numpy as np

np.random.seed(2)

A = np.random.uniform(0, 10, (10, 2))

triangle_enclosing(A, L = 40)
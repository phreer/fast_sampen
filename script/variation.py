import math
from functools import cmp_to_key
import sys
import numpy as np

from utils import gen_points
from utils import distance_within

n = 1000
m = 3
points, var = gen_points(sys.argv[1], n, m)

def less_point(p1, p2):
    for i in range(m):
        if p1[i] < p2[i]: return -1
        elif p1[i] > p2[i]: return 1
    return 0

r = 0.1 * math.sqrt(var)

zuv = np.zeros([n, n], np.int32)
for i in range(len(points)):
    for j in range(i + 1, len(points)):
        if distance_within(points[i], points[j], r, m):
            zuv[i, j] = 1
            zuv[j, i] = 1

As = np.sum(zuv, axis=1)
vas = np.abs(As[1: ] - As[: -1]).sum()
vzuv = np.abs(zuv[1 : , : ] - zuv[ : -1, : ]).sum()

print('Result (no presort)')
print('A: ', As.sum())
print('Variation of {A_i}: ', vas)
print('Variation of {z_uv}: ', vzuv)

points = sorted(points, key=cmp_to_key(less_point))
zuv = np.zeros([n, n], np.int32)
for i in range(len(points)):
    for j in range(i + 1, len(points)):
        if distance_within(points[i], points[j], r, m):
            zuv[i, j] = 1
            zuv[j, i] = 1

As = np.sum(zuv, axis=1)
vas = np.abs(As[1: ] - As[: -1]).sum()
vzuv = np.abs(zuv[1 : , : ] - zuv[ : -1, : ]).sum()

print('Result (presort)')
print('A: ', As.sum())
print('Variation of {A_i}: ', vas)
print('Variation of {z_uv}: ', vzuv)

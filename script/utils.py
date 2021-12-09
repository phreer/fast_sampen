import numpy as np


def gen_points(filename, n, m):
    xs = []
    points = []
    with open(filename) as f:
        xs = [next(f).split()[1] for i in range(n + m)]
        xs = list(map(float, xs))
    xs = np.array(xs)
    var = np.sum((xs - np.mean(xs)) ** 2) / (len(xs) - 1)
    for i in range(n):
        p = xs[i: i + m]
        points.append(p)
    return points, var


def distance_within(p1, p2, r, m):
    for i in range(m):
        if p1[i] - p2[i] > r or p1[i] - p2[i] < -r:
            return False
    return True

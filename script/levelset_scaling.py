import numpy as np
import matplotlib.pyplot as plt
from freudenthal import freudenthal_grid
import bats
import gudhi as gd
#import dionysus as d
import time

def time_BATS_update(img, img2):
    m, n = img.shape
    X = freudenthal_grid(m, n)
    
    t0 = time.monotonic()
    vals, imap = bats.lower_star_filtration(X, img.flatten())
    t1 = time.monotonic()
    print("time to extend: {} sec.".format(t1 - t0))
    
    t0 = time.monotonic()
    F = bats.FilteredSimplicialComplex(X, vals)
    t1 = time.monotonic()
    print("time to construct: {} sec.".format(t1 - t0))
    
    t0 = time.monotonic()
    R = bats.reduce(F, bats.F2())
    t1 = time.monotonic()
    print("time to reduce: {} sec.".format(t1 - t0))
    
    t0 = time.monotonic()
    vals, imap = bats.lower_star_filtration(X, img2.flatten())
    R.update_filtration(vals)
    t1 = time.monotonic()
    print("time to update: {} sec.".format(t1 - t0))
    tupd = t1 - t0
    
    t0 = time.monotonic()
    vals, imap = bats.lower_star_filtration(X, img2.flatten())
    F = bats.FilteredSimplicialComplex(X, vals)
    R = bats.reduce(F, bats.F2())
    t1 = time.monotonic()
    print("img2 from scratch: {} sec.".format(t1 - t0))
    tscratch = t1 - t0

    return tupd, tscratch


tupd = []
tscratch = []
ns = [50, 100, 150, 200, 250, 300, 350, 400]
for n in ns:
    print("\n{}".format(n))
    img = np.empty((n,n), dtype=np.float64)
    for i in range(n):
        for j in range(n):
            img[i,j] = np.sin(10* np.pi * i / n) + np.cos(10* np.pi * j/n)
    
    img2 = img + 0.01 * np.random.randn(n,n)
    
    ta, tb = time_BATS_update(img, img2)
    tupd.append(ta)
    tscratch.append(tb)


print('n')
print(ns)
print('tupd')
print(tupd)
print('tscratch')
print(tscratch)

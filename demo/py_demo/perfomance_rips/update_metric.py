import bats
import numpy as np
import time

n = 100
X = np.random.normal(size=(n,2))

def time_BATS_update_L1(X):
    tlist = []
    
    m1 = bats.Euclidean()
    m2 = bats.L1Dist()
    
    data = bats.DataSet(bats.Matrix(X)) # put into a bats.DataSet
    
    t0 = time.monotonic()
    F_X = bats.RipsFiltration(data, m1, np.inf, 2) # generate a RipsFiltration
    R_X = bats.reduce(F_X, bats.F2()) # reduce with F2 coefficients
    t1 = time.monotonic()
    tlist.append(t1-t0)
    #print("time to reduce: {} sec.".format(t1 - t0))
    
    t0 = time.monotonic()
    F_Y = bats.RipsFiltration(data, m2, np.inf, 2) # generate a new RipsFiltration
    UI = bats.UpdateInfoFiltration(F_X, F_Y)
    R_X.update_filtration_general(UI)
    t1 = time.monotonic()
    #print("time to update: {} sec.".format(t1 - t0))
    tlist.append(t1-t0)
    
    return tlist


record = []
for i in range(30):
    record.append(time_BATS_update_L1(X))

# get mean time
record = np.array(record)
np.mean(record, axis = 0)

import gudhi as gd
t0 = time.monotonic()
rips_complex = gd.RipsComplex(points=X, max_edge_length=np.inf)
simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
diag = simplex_tree.persistence()
t1 = time.monotonic()
print("compute1: {} sec.".format(t1 - t0))

t0 = time.monotonic()
m1 = bats.Euclidean()
data = bats.DataSet(bats.Matrix(X)) # put into a bats.DataSet
F_X = bats.RipsFiltration(data, m1, np.inf, 2) # generate a RipsFiltration
R_X = bats.reduce(F_X, bats.F2()) # reduce with F2 coefficients
t1 = time.monotonic()
print("compute1: {} sec.".format(t1 - t0))
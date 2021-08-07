import bats
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.spatial.distance as distance
import csv

import torch
from topologylayer.nn import RipsLayer # topologylayer package
import gudhi as gd # Gudhi
import dionysus
from ripser import ripser

# Don't run this for larger problems
def time_toplayer_rips(y, dmax = 1):

    t0 = time.monotonic()
    layer = RipsLayer(y.shape[0], maxdim = dmax) # maxdim is homology dimension
    y_t = torch.tensor(y, dtype=torch.float, requires_grad=False)
    pds = layer(y_t)
    t1 = time.monotonic()

    return t1 - t0

def time_gudhi_rips(X, dmax):

    t0 = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    rips_complex = gd.RipsComplex(points=X, max_edge_length=rX)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension = dmax)
    diag = simplex_tree.persistence()
    t1 = time.monotonic()

    return t1 - t0

def time_dionysus_rips(X, dmax):

    t0 = time.monotonic()
    # construct
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    f = dionysus.fill_rips(X, dmax, rX)
    p = dionysus.homology_persistence(f)  # using clearing method by default
    t1 = time.monotonic()

    return t1 - t0

def time_ripser_rips(X, **kw):

    t0 = time.monotonic()
    out = ripser(X, **kw)
    t1 = time.monotonic()
    return t1 - t0


# time to compute rips in BATs with flags,
# Return total time
def time_BATS_rips_flags(X, dmax = 2, flags=(bats.standard_reduction_flag(), bats.compute_basis_flag())):

    t0 = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    F = bats.LightRipsFiltration(bats.Matrix(DX), rX, dmax)
    R = bats.reduce(F, bats.F2(), *flags)
    t1 = time.monotonic()
    #print("reduce1: {} sec.".format(t1 - t0))

    return t1-t0



# update persistence on Y and only return total time on updating
def time_BATS_updates_rips(X, Y, dmax = 2):

    # compution on X
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    F = bats.LightRipsFiltration(bats.Matrix(DX), rX, dmax)
    R = bats.reduce(F, bats.F2())

    # return time spent on updating Y
    t0 = time.monotonic()
    DY = distance.squareform(distance.pdist(Y))
    rY = bats.enclosing_radius(bats.Matrix(DY))
    FY = bats.LightRipsFiltration(bats.Matrix(DY), rY, dmax)
    # update
    update_info = bats.UpdateInfoLightFiltration(F, FY)
    R.update_filtration_general(update_info)
    t1 = time.monotonic()

    return t1-t0



# generate circle (unnoisy)
def gen_circle(n = 50, d = 2):
    X = np.random.normal(size=(n,d))
    # X = np.random.uniform(-1,1,(n,2))
    X = X / np.linalg.norm(X, axis=1).reshape(-1,1)
    # std = 0.1
    # Y = Y + np.random.normal(size=(n,2), scale = 0.1 )
    return X


# Now ready to run
# write into files
header = ['dataset',
'standard_basis',
'standard_no_basis',
'clearing',
'compression',
'update',
'topologylayer',
'gudhi',
'dionysus',
'ripser'] # set header of csv file
data = [] # create a list used to write into file



flags = [
    (bats.standard_reduction_flag(), bats.compute_basis_flag()),
    (bats.standard_reduction_flag(),),
    (bats.standard_reduction_flag(), bats.clearing_flag()),
    (bats.standard_reduction_flag(), bats.compression_flag())
]
labels = [
    "standard w/ basis",
    "standard w/ no basis",
    "standard w/ clearing",
    "standard w/ compression",
]




# check if header is needer to add into file
need_header = True
import os.path
filename = 'bats_comp_other_package.csv'
if os.path.isfile(filename):
    need_header = False

with open(filename, 'a', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    if need_header:
        writer.writerow(header)

    # noisy circle dataset
    for sigma in [1e-3]:
        for n in [20]:
            for d in [3]:
                # create dataset
                X = gen_circle(n, d) # unnoisy circle
                # X = X + np.random.normal(size=(n,2), scale = 0.1)
                # Y = X + 0.005*np.random.randn(n,2)
                X2 = X + np.random.normal(size=(n,d), scale = sigma )

                # for each dataset repeat
                n_repeats = 1
                for i in range(n_repeats):
                    # write row
                    record_line = [] # record for each line
                    s1 = 'sphere-'
                    s2 = str(sigma) + '-'
                    s3 = str(n) + '-'
                    s4 = str(d)
                    s5 = s1 + s2 + s3 + s4
                    record_line.append(s5) # append dataset name

                    for flag, label in zip(flags, labels):
                        record_line.append(time_BATS_rips_flags(X2, d, flag))

                    record_line.append(time_BATS_updates_rips(X, X2, d))
                    if n > 30 or d > 2: # topologylayer is slow and so we do not run on larger problems
                        record_line.append(0)
                    else:
                        record_line.append(time_toplayer_rips(X2,d))
                    record_line.append(time_gudhi_rips(X2, d))
                    record_line.append(time_dionysus_rips(X2,d))
                    record_line.append(time_ripser_rips(X2,maxdim = d))
                    # add row into data
                    data.append(record_line)

    # print("with basis")
    # time_ripser_rips(Y, Y2, do_cocycles=True, maxdim=1)


    # write multiple rows
    writer.writerows(data)

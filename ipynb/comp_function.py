# this file will define functions used to compare time between packages
# also a function to return a record line in csv file
import numpy as np
import pandas as pd
import time
import scipy.spatial.distance as distance
import bats
import gudhi as gd # Gudhi
import dionysus
from ripser import ripser

# import torch
# from topologylayer.nn import RipsLayer # topologylayer package
# # Don't run this for larger problems
# # Toplayer from
# # https://github.com/bruel-gabrielsson/TopologyLayer#ripslayer
# def time_toplayer_rips(y, Hmaxd = 1):

#     t0 = time.monotonic()
#     layer = RipsLayer(y.shape[0], maxdim = Hmaxd) # maxdim is homology dimension
#     y_t = torch.tensor(y, dtype=torch.float, requires_grad=False)
#     pds = layer(y_t)
#     t1 = time.monotonic()

#     return t1 - t0

# GUDHI from
# https://gudhi.inria.fr/python/latest/rips_complex_user.html
def time_gudhi_rips(X, dmax):
    # dmax is the dimension of complex
    t0 = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    rips_complex = gd.RipsComplex(points=X, max_edge_length=rX)
    simplex_tree = rips_complex.create_simplex_tree(max_dimension = dmax)
    diag = simplex_tree.persistence()
    t1 = time.monotonic()

    return t1 - t0

# Dionysus2 from
# https://mrzv.org/software/dionysus2/tutorial/rips.html
def time_dionysus_rips(X, dmax, degree = -1):
    # dmax is the dimension of complex
    t0 = time.monotonic()
    # construct
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    f = dionysus.fill_rips(X, dmax, rX)
    if degree == -1:
        p = dionysus.homology_persistence(f)  # using clearing method by default
    elif degree == 1:
        p = dionysus.cohomology_persistence(f, prime=2)
    else:
        print("Degree must be +1 or -1 but you type in", degree)
    t1 = time.monotonic()

    return t1 - t0

# Ripser from
# https://ripser.scikit-tda.org/en/latest/reference/stubs/ripser.ripser.html#ripser.ripser
def time_ripser_rips(X, **kw):
    t0 = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    out = ripser(X, thresh=rX, **kw)
    t1 = time.monotonic()

    return t1 - t0


# enc radius
# dmax is dimension of complex
# will compue H_0, H_1, ..., until H_{dmax-1}
# default flag: clearing with basis
def time_BATS_rips_enc(X, dmax = 2, degree = +1,
flags = (bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())):
    t0 = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    F_X = bats.LightRipsFiltration(bats.Matrix(DX), rX, dmax)
    FVS = bats.FilteredF2DGVectorSpace(F_X, degree)
    RFVS = bats.ReducedFilteredF2DGVectorSpace(FVS, *flags)
    t1 = time.monotonic()
    return t1-t0


# Update with max radius (see the effect of permutations)
# Note: update must include flag bats.compute_basis_flag()
# dmax is max dimension of complex
def time_BATS_update_full_rips(X, Y, rX = np.inf, rY = np.inf, dmax = 2, degree = +1,
flags = (bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())):
    # computation on X
    F_X = bats.LightRipsFiltration(bats.DataSet(bats.Matrix(X)), bats.Euclidean(), rX, dmax)
    FVS = bats.FilteredF2DGVectorSpace(F_X, degree)
    RFVS = bats.ReducedFilteredF2DGVectorSpace(FVS, *flags) 

    # return time spent on updating Y
    t0 = time.monotonic()
    F_Y = bats.LightRipsFiltration(bats.DataSet(bats.Matrix(Y)), bats.Euclidean(), rY, dmax)
    UI = bats.UpdateInfo2(F_X, F_Y)
    RFVS.update(UI)
    t1 = time.monotonic()
    # print("Update: {} sec.".format(t1 - t0))
    return t1-t0

# Update with enc radius (see the effect of permutations/addition/deletion)
def time_BATS_update_enc_rips(X, Y, dmax = 2, degree = +1,
flags = (bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())):
    # computation on X
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    F_X = bats.LightRipsFiltration(bats.Matrix(DX), rX, dmax)
    FVS = bats.FilteredF2DGVectorSpace(F_X, degree)
    RFVS = bats.ReducedFilteredF2DGVectorSpace(FVS, *flags) 

    # return time spent on updating Y
    t0 = time.monotonic()
    DY = distance.squareform(distance.pdist(Y))
    rY = bats.enclosing_radius(bats.Matrix(DY))
    F_Y = bats.LightRipsFiltration(bats.Matrix(DY), rY, dmax)
    UI = bats.UpdateInfo2(F_X, F_Y)
    RFVS.update(UI)
    t1 = time.monotonic()
    # print("Update: {} sec.".format(t1 - t0))
    return t1-t0

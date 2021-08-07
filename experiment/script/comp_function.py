# this file will define functions used to compare time between packages
# also a function to return a record line in csv file
from plyfile import PlyData, PlyElement
import numpy as np
import pandas as pd
import time
import scipy.spatial.distance as distance


import bats
import torch
from topologylayer.nn import RipsLayer # topologylayer package
import gudhi as gd # Gudhi
import dionysus
from ripser import ripser


# Don't run this for larger problems
# Toplayer from
# https://github.com/bruel-gabrielsson/TopologyLayer#ripslayer
def time_toplayer_rips(y, Hmaxd = 1):

    t0 = time.monotonic()
    layer = RipsLayer(y.shape[0], maxdim = Hmaxd) # maxdim is homology dimension
    y_t = torch.tensor(y, dtype=torch.float, requires_grad=False)
    pds = layer(y_t)
    t1 = time.monotonic()

    return t1 - t0

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
def time_dionysus_rips(X, dmax):
    # dmax is the dimension of complex
    t0 = time.monotonic()
    # construct
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    f = dionysus.fill_rips(X, dmax, rX)
    p = dionysus.homology_persistence(f)  # using clearing method by default
    t1 = time.monotonic()

    return t1 - t0

# Ripser from
# https://ripser.scikit-tda.org/en/latest/reference/stubs/ripser.ripser.html#ripser.ripser
def time_ripser_rips(X, **kw):
    t0 = time.monotonic()
    out = ripser(X, **kw)
    t1 = time.monotonic()

    return t1 - t0


# time to compute rips in BATs with flags,
# Return total time
def time_BATS_rips_flags(X, dmax = 2, flags=(bats.standard_reduction_flag(), bats.compute_basis_flag())):
    # dmax is dimension of complex
    # will compue H_0, H_1, ..., until H_{dmax-1}
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


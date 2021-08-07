import bats
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.spatial.distance as distance
import csv


# time to compute rips in BATs with different flags,
# Return total time
def time_BATS_rips_flags(X, dmax = 2, flags=(bats.standard_reduction_flag(), bats.compute_basis_flag())):
    t0 = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    F = bats.LightRipsFiltration(bats.Matrix(DX), rX, dmax)
    R = bats.reduce(F, bats.F2(), *flags)
    t1 = time.monotonic()

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
header = ['std','n_points', 'dimension',
'standard_basis',
'standard_no_basis',
'clearing',
'compression',
'update'] # set header of csv file
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
filename = 'bats_flag_compare_to_updating.csv'
if os.path.isfile(filename):
    need_header = False

with open(filename, 'a', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    if need_header:
        writer.writerow(header)

    for sigma in [2e-3]:
        for n in [100]:
            for d in [2,3]:
                # create dataset
                X = gen_circle(n, d) # unnoisy circle
                # X = X + np.random.normal(size=(n,2), scale = 0.1)
                # Y = X + 0.005*np.random.randn(n,2)
                X2 = X + np.random.normal(size=(n,d), scale = sigma )

                # for each dataset repeat
                n_repeats = 20
                for i in range(n_repeats):
                    # write row
                    record_line = [] # record for each line
                    record_line.append(sigma)
                    record_line.append(n)
                    record_line.append(d)

                    for flag, label in zip(flags, labels):
                        record_line.append(time_BATS_rips_flags(X2, d, flag))

                    record_line.append(time_BATS_updates_rips(X, X2, d))

                    # add row into data
                    data.append(record_line)

    # write multiple rows
    writer.writerows(data)

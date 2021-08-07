# Generate a noisy circle and then add normal noise with std = 0.005 to it
# Experiment on different size of datasets

import bats
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.spatial.distance as distance
import csv

# first give a function to run

def time_BATS_updates_enc_rips_record_from_euc_to_minkow1(X, dmax = 2):
    rec = []

    t0 = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    t1 = time.monotonic()
    #print("setup1: {} sec.".format(t1 - t0))
    rec.append(t1-t0)

    t0 = time.monotonic()
    F = bats.LightRipsFiltration(bats.Matrix(DX), rX, dmax)
    t1 = time.monotonic()
    #print("construct1: {} sec.".format(t1 - t0))
    rec.append(t1-t0)

    t0 = time.monotonic()
    R = bats.reduce(F, bats.F2())
    t1 = time.monotonic()
    #print("reduce1: {} sec.".format(t1 - t0))
    rec.append(t1-t0)
    nnz_old_U = np.sum(R.nnz_U())
    nnz_old_R = np.sum(R.nnz_R())

    t0 = time.monotonic()
    DY = distance.squareform(distance.pdist(X, 'minkowski', p=1))
    rY = bats.enclosing_radius(bats.Matrix(DY))
    t1 = time.monotonic()
    #print("setup2: {} sec.".format(t1 - t0))
    rec.append(t1-t0)

    t0 = time.monotonic()
    FY = bats.LightRipsFiltration(bats.Matrix(DY), rY, dmax)
    t1 = time.monotonic()
    #print("construct2: {} sec.".format(t1 - t0))
    rec.append(t1-t0)

    t0 = time.monotonic()
    RY = bats.reduce(FY, bats.F2())
    t1 = time.monotonic()
    #print("reduce2: {} sec.".format(t1 - t0))
    rec.append(t1-t0)
    # t_reb = t1-t0
    nnz_reb_U = np.sum(RY.nnz_U())
    nnz_reb_R = np.sum(RY.nnz_R())

    t0 = time.monotonic()
    update_info = bats.UpdateInfoLightFiltration(F, FY)
    R.update_filtration_general(update_info)
    t1 = time.monotonic()
    #print("update: {} sec.".format(t1 - t0))
    rec.append(t1-t0)
    # t_upd = t1-t0
    # factor = round(t_reb/t_upd, 4)
    # rec.append(factor)
    nnz_upd_U = np.sum(R.nnz_U())
    nnz_upd_R = np.sum(R.nnz_R())

    rec.append(nnz_old_U)
    rec.append(nnz_reb_U)
    rec.append(nnz_upd_U)
    rec.append(nnz_old_R)
    rec.append(nnz_reb_R)
    rec.append(nnz_upd_R)

    return rec

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
'setup1', 'construct1', 'reduce1',
'setup2', 'construct2', 'reduce2',
'update',
'nnz_old_U', 'nnz_reb_U', 'nnz_upd_U',
'nnz_old_R', 'nnz_reb_R', 'nnz_upd_R'] # set header of csv file
data = [] # create a list used to write into file

# check if header is needer to add into file
need_header = True
import os.path
filename = 'noisy_circle_update_metric_euc_to_minkow1.csv'
if os.path.isfile(filename):
    need_header = False

with open(filename, 'a', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    if need_header:
        writer.writerow(header)

    for sigma in [0.005]:
        for n in [100]:
            for d in [2]:
                # create dataset
                X = gen_circle(n, d) # unnoisy circle
                # X = X + np.random.normal(size=(n,2), scale = 0.1)
                # Y = X + 0.005*np.random.randn(n,2)
                X = X + np.random.normal(size=(n,d), scale = sigma )

                # for each dataset repeat
                n_repeats = 20
                for i in range(n_repeats):
                    record_line = [] # record for each line
                    record_line.append(sigma)
                    record_line.append(n)
                    record_line.append(d)
                    record_line.extend(time_BATS_updates_enc_rips_record_from_euc_to_minkow1(X,d))

                    # add row into data
                    data.append(record_line) # add row

    # write multiple rows
    writer.writerows(data)

"""
This is a file that records the time used for updating Rips Filtration 
from Euclidean to L1 metric.
"""


import bats
import numpy as np
import time # used to compute time
import csv # used to store csv file

# update metric from Euclidean to L1 
def time_BATS_update_L1(X):
    tlist = []
    
    m1 = bats.Euclidean()
    m2 = bats.L1Dist()
    dmax = 3
    rmax = np.inf
    data = bats.DataSet(bats.Matrix(X)) # put into a bats.DataSet
    
    t0 = time.monotonic()
    F_X = bats.RipsFiltration(data, m1, rmax, dmax) # generate a RipsFiltration
    R_X = bats.reduce(F_X, bats.F2()) # reduce with F2 coefficients
    t1 = time.monotonic()
    tlist.append(t1-t0)
    #print("time to reduce: {} sec.".format(t1 - t0))
    
    # Rebuild
    t0 = time.monotonic()
    F_Y = bats.RipsFiltration(data, m2, rmax, dmax) # generate a RipsFiltration
    R_Y = bats.reduce(F_X, bats.F2()) # reduce with F2 coefficients
    t1 = time.monotonic()
    tlist.append(t1-t0)

    # Update
    t0 = time.monotonic()
    F_Y = bats.RipsFiltration(data, m2, rmax, dmax) # generate a new RipsFiltration
    UI = bats.UpdateInfoFiltration(F_X, F_Y)
    R_X.update_filtration_general(UI)
    t1 = time.monotonic()
    #print("time to update: {} sec.".format(t1 - t0))
    tlist.append(t1-t0)
    
    return tlist




# write into files
header = ['initial_redue', 'rebuild', 'update'] # set header of csv file
data = []

with open('euc_to_L1.csv', 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    for n in [20, 40]:

        # create dataset
        X = np.random.normal(size=(n,2))
        record_n = []
        num_repeats = 2
        for i in range(num_repeats):
            record_n.append(time_BATS_update_L1(X))

        # get mean time and add to records
        record_n = np.array(record_n)
        data.append(np.mean(record_n, axis = 0))

    # write multiple rows
    writer.writerows(data)


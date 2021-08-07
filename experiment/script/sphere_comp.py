import os
from plyfile import PlyData, PlyElement
import numpy as np
import pandas as pd
import bats
import matplotlib.pyplot as plt
import time
import scipy.spatial.distance as distance
import csv
# load functions that will record time
from comp_function import time_BATS_rips_flags
from comp_function import time_BATS_updates_rips
from comp_function import time_dionysus_rips
from comp_function import time_gudhi_rips
from comp_function import time_ripser_rips
from comp_function import time_toplayer_rips

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
header = ['n_points','dimension','std','standard_basis',
'standard_no_basis',
'clearing',
'compression',
'update',
'topologylayer',
'gudhi',
'dionysus',
'ripser'] # set header of csv file



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

filename = 'sphere_record100_H1.csv'
with open(filename, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    data_into_file = [] # create a list used to write into file
    # noisy circle dataset
    for sigma in [1e-3]:
        for n in [100]:
            for d in [2]:
                # set parameter
                n_repeats = 20
                dmax_complex = d # maximum dimension of complex, will compute PH until one dimension lower
                for i in range(n_repeats):
                    # create dataset
                    X = gen_circle(n, d) # unnoisy circle
                    # X = X + np.random.normal(size=(n,2), scale = 0.1)
                    # Y = X + 0.005*np.random.randn(n,2)
                    X2 = X + np.random.normal(size=(n,d), scale = sigma )

                    # write row
                    record_line = [] # record for each line
                    record_line.append(n)
                    record_line.append(d)
                    record_line.append(sigma)

                    # BATs Method with Flags
                    for flag, label in zip(flags, labels):
                        record_line.append(time_BATS_rips_flags(X2, dmax_complex, flag))

                    # BATs update method
                    record_line.append(time_BATS_updates_rips(X, X2, dmax_complex))

                    # topologylayer
                    record_line.append(0) # topologylayer is slow not run
                    # record_line.append(time_toplayer_rips(X2,dmax_complex-1))

                    # Gudhi
                    record_line.append(time_gudhi_rips(X2, dmax_complex))
                    
                    # dionysus
                    # record_line.append(0) # dionysus slow not run
                    record_line.append(time_dionysus_rips(X2, dmax_complex)) 

                    # ripser
                    # Take care, maxdim is Maximum homology dimension
                    # dmax_complex is the complex dimension, need to mius one
                    record_line.append(time_ripser_rips(X2, maxdim = dmax_complex-1)) 

                    # add row into data
                    data_into_file.append(record_line)



    # write multiple rows
    writer.writerows(data_into_file)

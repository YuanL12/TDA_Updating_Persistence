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





# Now ready to run
# write into files
header = ['standard_basis','standard_no_basis','clearing','compression','update',
'topologylayer','gudhi','dionysus','ripser'] # set header of csv file



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


# csv file name
dmax_complex = 3 # maximum dimension of complex, will compute PH until one dimension lower
n_sample = 100 # number of sample points
filename = 'klein_record{}_H{}.csv'.format(n_sample, dmax_complex-1)

with open(filename, 'w', encoding='UTF8', newline='') as f:
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    # create dataset
    df = pd.read_csv('data/roadmap_datasets_point_cloud/klein_bottle_pointcloud_new_400.txt', delimiter = " ", header = None)
    df = df.iloc[:,:-1] # the last column is NaN, we need it
    
    
    # set parameter
    n_repeats = 20
    data_into_file = [] # create a list used to write into file
    for i in range(n_repeats):

        X = np.array(df) # transform into numpy array
        # Uniform randomly generate points
        n_total_rows = X.shape[0]
        random_indices = np.random.choice(n_total_rows, size=n_sample, replace=False)
        X = X[random_indices, :]
        # X2 = X + X * 0.01 * np.random.randn(X.shape[0],X.shape[1])
        X2 = X + np.random.normal(size=(X.shape[0], X.shape[1]), scale = 0.001 )

        # write row
        record_line = [] # record for each line

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
        record_line.append(0) # dionysus slow not run
        # record_line.append(time_dionysus_rips(X2, d)) # dionysus slow not run
        
        # ripser
        # Take care, maxdim is Maximum homology dimension
        # dmax_complex is the complex dimension, need to mius one
        record_line.append(time_ripser_rips(X2, maxdim = dmax_complex-1)) 

        # add row into data
        data_into_file.append(record_line)


    # write multiple rows
    writer.writerows(data_into_file)

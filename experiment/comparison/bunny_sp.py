# Subsamplying for bunny dataset
import numpy as np
import bats
from bats.visualization.plotly import MapVisualization
import scipy.spatial.distance as distance
import os
from plyfile import PlyData, PlyElement
import pandas as pd


def read_ply_to_numpy_arr(file_path):

    plydata = PlyData.read(file_path) # read file
    data = plydata.elements[0].data # read data
    data_pd = pd.DataFrame(data) # Convert to DataFrame, because DataFrame can parse structured data
    data_np = np.zeros(data_pd.shape, dtype=np.float32) # Initialize the array of stored data
    property_names = data[0].dtype.names # read the name of the property
    for i, name in enumerate(property_names): # Read data according to property, so as to ensure that the data read is the same data type.
        data_np[:, i] = data_pd[name]

    return data_np[:,:3] # only use the first 3 columns

X = read_ply_to_numpy_arr('bunny/data/bun000.ply')

DX = distance.squareform(distance.pdist(X, 'euclidean'))
# rX = bats.enclosing_radius(bats.Matrix(DX))
# R = bats.RipsComplex(bats.Matrix(DX), rX, 2)

# store distance matrix
pd.DataFrame(DX).to_csv("pd_bunny.csv",header = None, index = None)
# df = pd.read_csv("pd_cylider.csv", header = None)

# inds is sequence of indices selected by greedy landmarking
# dists[k] is hausdorff distance from X[inds[:k]] to X
inds, dists = bats.greedy_landmarks_hausdorff(bats.Matrix(DX),0) # 0 is first index

for k in [100,200,500,1000]:
    # we'll landmark k points
    # print('hausdorff distance is {}'.format(dists[k]))
    # ret = plt.scatter(X[inds[:k],0], X[inds[:k],1])
    X_sp = X[inds[:k],:]
    # use scipy to get pairwise distances
    DX_sp = distance.squareform(distance.pdist(X_sp, 'euclidean'))
    pd.DataFrame(DX_sp).to_csv("pd_bunny_sp{}_HD_{}.csv".format(k, dists[k]), header = None, index = None)
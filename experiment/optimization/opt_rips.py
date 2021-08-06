# OPTIMIZE: increase \epsilon (2,0,1;PD1)


import csv
import bats
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import scipy.spatial.distance as distance


# time to compute rips in BATs with flags,
# Return total time
def time_BATS_rips_flags(X, dmax = 2, flags=(bats.standard_reduction_flag(), bats.compute_basis_flag())):
    t0 = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    FCC, ip = bats.LightRipsFiltration_extension(bats.Matrix(DX), rX, dmax)
    RFCC = bats.reduce(FCC, bats.F2(), *flags)
    t1 = time.monotonic()
    #print("reduce1: {} sec.".format(t1 - t0))

    return t1-t0, RFCC



n = 100
X = np.random.uniform(0,1,(n,2))
pd.DataFrame(X).to_csv("../csv_record/original_dataset_rips.csv",header = None, index = None)
# X = X / np.linalg.norm(X, axis=1).reshape(-1,1)
# X = X + np.random.normal(size=(n,2), scale = 0.1 )
# fig1 = plt.scatter(X[:,0], X[:,1])
# fig1.axes.set_aspect('equal')
# plt.savefig('original_dataset_rips.svg')
# plt.clf()


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

nnz_ = []
tlist = []
lr = 1e-2 # persistence penalty


t0 = time.monotonic()
DX = distance.squareform(distance.pdist(X))
rX = bats.enclosing_radius(bats.Matrix(DX))
F, imap = bats.LightRipsFiltration_extension(bats.Matrix(DX), rX , 2)
R = bats.reduce(F, bats.F2())
t1 = time.monotonic()
#print("initialization: {} sec.".format(t1-t0))
F_old = F

iter_times = 100
for i in range(iter_times):
    t_iter = []
    nnz_iter = []
    # print("iter", i)
    # Optimize H1 length
    # get 1-dimensional pairs
    t0 = time.monotonic()
    ps = R.persistence_pairs(1)
    t1a = time.monotonic()
    #print("\tfind pairs: {} sec.".format(t1a - t0))

    t0b = time.monotonic()
    complex_old = F_old.complex() # complex of filtration
    filtration_vals_old = F_old.vals(1) # filtration vals at dimension 1

    # change the posistion of two vertices
    # that are related to the birth of the persistence pair

    #print("rX = ", rX)
    for p in ps:
        if p.length() > 0:
        # if p.length() > 0 and np.abs(p.birth() - rX) > 1e-4:
            #print(p)
            d = p.dim()
            # get indices of two vertices of the (birth) edge
            bi = p.birth_ind() # index of birth edge
            [birth_vertex1_i, birth_vertex2_i] = complex_old.get_simplex(1, bi)
            # Gradient Descent
            grad_div_2 = 2 * (p.death() - p.birth())* (X[birth_vertex1_i] - X[birth_vertex2_i])
            X[birth_vertex1_i] -=  (lr/ filtration_vals_old[bi]) * (grad_div_2)
            X[birth_vertex2_i] +=  (lr/ filtration_vals_old[bi]) * (grad_div_2)
            #print("update birth success")
            # get the death index of related edge
            if p.death_ind() > len(imap[d+1])-1:
                print("death index exceed limit: ", p.death_ind())
                print("p is ", p)

            if p.death() != float('inf') and p.death_ind() <= len(imap[d+1])-1:
                di = imap[d+1][p.death_ind()] # maps death_ind to the death edge (related to the 2-simplex destroys H1)
                # get index of two vertices of the (birth) edge
                [death_vertex1_i, death_vertex2_i] = complex_old.get_simplex(1, di)
                # Gradient Descent
                grad_div_2 = 2 * (p.death() - p.birth()) *  (X[death_vertex1_i] - X[death_vertex2_i])
                X[death_vertex1_i] +=  (lr/ filtration_vals_old[di]) * (grad_div_2)
                X[death_vertex2_i] -=  (lr/ filtration_vals_old[di]) * (grad_div_2)
                #print("update death success")


    t1b = time.monotonic()
    #print("\tUpdate Dataset: {} sec.".format(t1b - t0b))


    # Rebuild with different flags
    # Notice that only with basis can return the number of non-zeros
    # in U!!!
    for flag in flags:
        t_reb, R_reb = time_BATS_rips_flags(X, 2, flag)
        t_iter.append(t_reb)
        # nnz_iter.append(np.sum(R_reb.nnz_U() + R_reb.nnz_R())) will return error
        nnz_iter.append(np.sum(R_reb.nnz_R()))


    # Update
    t0c = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    F_new, imap = bats.LightRipsFiltration_extension(bats.Matrix(DX), rX, 2)
    UI = bats.UpdateInfoLightFiltration(F_old, F_new)
    R.update_filtration_general(UI)
    t1c = time.monotonic()

    # Record updating time and nnz
    nnz_iter.append(np.sum(R.nnz_U() + R.nnz_R()))
    t_iter.append(t1c - t0c)

    tlist.append(t_iter)
    nnz_.append(nnz_iter)

    # Set old Filtration
    F_old = F_new



# nnz are only for R matrix!!!!!
header = ["std_basis", "std_no_basis", "std_clearing", "std_compression",'update',
'nnz_basis','nnz_no_basis', 'nnz_clear','nnz_compress','nnz_upd'] # set header of csv file

data = []


# check if header is needer to add into file
# need_header = True
# import os.path
filename = '../csv_record/opt_enc_rips.csv'
# if os.path.isfile(filename):
#     need_header = False

with open(filename, 'w', encoding='UTF8', newline='') as f: # always rewrite the csv file
    writer = csv.writer(f)

    # write the header
    writer.writerow(header)

    data = np.concatenate((tlist, nnz_), axis=1)

    # write multiple rows
    writer.writerows(data)


# save optimized dataset
pd.DataFrame(X).to_csv("../csv_record/optimized_dataset_rips.csv", header = None, index = None)
# fig2 = plt.scatter(X[:,0], X[:,1])
# fig2.axes.set_aspect('equal')
# plt.savefig('optimized_dataset_rips.svg')
# plt.clf()

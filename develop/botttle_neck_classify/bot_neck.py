import bats
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.spatial.distance as distance
import gudhi as gd

n = 50
total_X = []
for rad in range(1,10):
    X = np.random.uniform(-1,1,(n,2))
    X = (rad* X) / np.linalg.norm(X, axis=1).reshape(-1,1)
    X = X + np.random.normal(size=(n,2), scale = 0.1 )
    total_X.append(X)
    fig1 = plt.scatter(X[:,0], X[:,1])
plt.savefig('orginal_dataset_rips.png')
plt.clf()

nnz_ = []
tlist = []
lr = 1e-2 # persistence penalty

t0 = time.monotonic()
X = total_X[0]
DX = distance.squareform(distance.pdist(X))
rX = bats.enclosing_radius(bats.Matrix(DX))
F, imap = bats.LightRipsFiltration_extension(bats.Matrix(DX), rX, 2)
R = bats.reduce(F, bats.F2())
t1 = time.monotonic()
print("Using enc_rad: {} sec.".format(t1-t0))

ps = R.persistence_pairs(0, False) +  R.persistence_pairs(1, False)
# bats.persistence_diagram(ps)
# plt.savefig('pd.png')
# plt.clf()
# print("enc_rad is ", rX)

total_nzps = []
nzps = []
for p in ps:
    if p.length() > 0:
        #print(p)
        bd_pair = []
        bd_pair.append(p.birth())
        bd_pair.append(p.death())
        nzps.append(bd_pair)

total_nzps.append(nzps)

for X in total_X[1:]:

    DX = distance.squareform(distance.pdist(X))
    rX = bats.enclosing_radius(bats.Matrix(DX))
    F = bats.LightRipsFiltration(bats.Matrix(DX), rX, 2)
    R = bats.reduce(F, bats.F2())

    ps = R.persistence_pairs(0, False) +  R.persistence_pairs(1, False)
    nzps = []
    for p in ps:
        if p.length() > 0:
            #print(p)
            bd_pair = []
            bd_pair.append(p.birth())
            bd_pair.append(p.death())
            nzps.append(bd_pair)

    total_nzps.append(nzps)

# Reconstruct a noisy circle with radius 4
X = np.random.uniform(-1,1,(n,2))
X = (4 * X) / np.linalg.norm(X, axis=1).reshape(-1,1)
X = X + np.random.normal(size=(n,2), scale = 0.1 )
DX = distance.squareform(distance.pdist(X))
rX = bats.enclosing_radius(bats.Matrix(DX))
F, imap = bats.LightRipsFiltration_extension(bats.Matrix(DX), rX, 2)
R = bats.reduce(F, bats.F2())

ps = R.persistence_pairs(0, False) +  R.persistence_pairs(1, False)

nzps = []
for p in ps:
    if p.length() > 0:
        #print(p)
        bd_pair = []
        bd_pair.append(p.birth())
        bd_pair.append(p.death())
        nzps.append(bd_pair)


# find the smallest distance
dists = []
for i in range(len(total_nzps)):
    dists.append(gd.bottleneck_distance(nzps, total_nzps[i]))

print("The closest radius is")
print(dists.index(min(dists)) + 1) 

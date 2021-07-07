import bats
import numpy as np
import matplotlib.pyplot as plt
import time
import scipy.spatial.distance as distance

n = 100
X = np.random.uniform(0,1,(100,2))
# X = X / np.linalg.norm(X, axis=1).reshape(-1,1)
# X = X + np.random.normal(size=(n,2), scale = 0.1 )
fig1 = plt.scatter(X[:,0], X[:,1])
plt.savefig('orginal_dataset_rips.png')
plt.clf()

nnz_ = []
tlist = []
lr = 1e-2 # persistence penalty


t0 = time.monotonic()
DX = distance.squareform(distance.pdist(X))
# rX = bats.enclosing_radius(bats.Matrix(DX))
F, imap = bats.LightRipsFiltration_extension(bats.Matrix(DX), np.inf , 2)
R = bats.reduce(F, bats.F2())
t1 = time.monotonic()
#print("initialization: {} sec.".format(t1-t0))
F_old = F

iter_times = 100
for i in range(iter_times):
    #print("iter", i)
    # Optimize H1 length
    # get 1-dimensional pairs
    t0 = time.monotonic()
    ps = R.persistence_pairs(1)
    t1a = time.monotonic()
    #print("\tpairs: {} sec.".format(t1a - t0))
    
    t0b = time.monotonic()
    complex_old = F_old.complex() # complex of filtration
    filtration_vals_old = F_old.vals(1) # filtration vals at dimension 1
    
    # change the posistion of two vertices 
    # that are related to the birth of the persistence pair
    for p in ps:
        if p.length() > 0:
            d = p.dim()
            # get indices of two vertices of the (birth) edge 
            bi = p.birth_ind() # index of birth edge 
            [birth_vertex1_i, birth_vertex2_i] = complex_old.get_simplex(1, bi)
            # Gradient Descent
            grad_div_2 = 2 * (p.death() - p.birth())* (X[birth_vertex1_i] - X[birth_vertex2_i])
            X[birth_vertex1_i] -=  (lr/ filtration_vals_old[bi]) * (grad_div_2)
            X[birth_vertex2_i] +=  (lr/ filtration_vals_old[bi]) * (grad_div_2)

            # get the death index of related edge 
            di = imap[d+1][p.death_ind()] # maps death_ind to the death edge (related to the 2-simplex destroys H1)
            # get index of two vertices of the (birth) edge 
            [death_vertex1_i, death_vertex2_i] = complex_old.get_simplex(1, di)
            # Gradient Descent
            grad_div_2 = 2 * (p.death() - p.birth())*  (X[death_vertex1_i] - X[death_vertex2_i])
            X[death_vertex1_i] +=  (lr/ filtration_vals_old[di]) * (grad_div_2)
            X[death_vertex2_i] -=  (lr/ filtration_vals_old[di]) * (grad_div_2)
        
    t1b = time.monotonic()
    #print("\tOptimize: {} sec.".format(t1b - t0b))
    
    # extend filtration
    t0c = time.monotonic()
    DX = distance.squareform(distance.pdist(X))
    F_new, imap = bats.LightRipsFiltration_extension(bats.Matrix(DX), np.inf, 2)
    t1c = time.monotonic()
    #print("\textension: {} sec.".format(t1c - t0c))
    
    # Find Updating Information
    t0d = time.monotonic()
    UI = bats.UpdateInfoLightFiltration(F_old, F_new)
    t1d = time.monotonic()
    #print("\tfind updating info: {} sec.".format(t1d - t0d))
    
    # Update persistence
    t0e = time.monotonic()
    R.update_filtration_general(UI)
    t1e = time.monotonic()
    #print("\tupdate: {} sec.".format(t1e - t0e))
    
    # Record updating time and nnz
    nnz_.append(np.sum(R.nnz_U() + R.nnz_R()))
    tlist.append(t1e - t0c)
    
    # Set old Filtration    
    F_old = F_new
    
    t1 = time.monotonic()
    #print("All: {} sec.".format(t1-t0))

fig2 = plt.scatter(X[:,0], X[:,1])
fig2.axes.set_aspect('equal')
plt.savefig('optimized_dataset_rips.png')
plt.clf()

fig3 = plt.figure(figsize=(8, 6))
plt.subplot(2, 1, 1)
plt.title('Sum of non-zeros in U and R v.s. number of iterations n')
plt.plot(np.arange(iter_times), nnz_, c = 'g')
plt.xlabel('number of iterations n')
plt.ylabel('number of non-zeros')


plt.subplot(2, 1, 2)
plt.title('Updating Time v.s. Sum of non-zeros of U and R')
plt.scatter(nnz_, tlist, c = 'r')
plt.xlabel('number of non-zeros')
plt.ylabel('Time')

fig3.tight_layout()
plt.savefig('opt_rips_results.png')
plt.clf()
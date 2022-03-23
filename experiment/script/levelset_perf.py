import bats
import gudhi as gd
import dionysus as dion
import ripser
import matplotlib.pyplot as plt
from freudenthal import *
import time
import numpy as np
from tqdm import tqdm

n = 128
img = np.empty((n,n), dtype=np.float64)
for i in range(n):
    for j in range(n):
        img[i,j] = np.sin(10* np.pi * i / n) + np.cos(10* np.pi * j/n)


m, n = img.shape
X = freudenthal_grid_light(m, n)

t0 = time.monotonic()
vals, imap = bats.lower_star_filtration(X, img.flatten())
t1 = time.monotonic()
print("time to extend: {} sec.".format(t1 - t0))

t0 = time.monotonic()
F = bats.FilteredLightSimplicialComplex(X, vals)
t1 = time.monotonic()
print("time to construct: {} sec.".format(t1 - t0))

t0 = time.monotonic()
R = bats.reduce(F, bats.F2())
t1 = time.monotonic()
print("time to reduce: {} sec.".format(t1 - t0))

# ps = R.persistence_pairs(0, False) +  R.persistence_pairs(1, False)
# for p in ps:
#     if p.length() > 0:
#         print(p)

N = 100
for sigma in [0.01, 0.1]:
    print("\n\n-----------------------------\n")
    print("variance: {}".format(sigma))
    data = sigma * np.random.randn(N, n, n) + img

    flags = [
        (),
        (bats.standard_reduction_flag(),),
        (bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag()),
        (bats.standard_reduction_flag(), bats.clearing_flag()),
        (bats.standard_reduction_flag(), bats.compression_flag()),
    ]
    labels = [
        "standard w/ basis",
        "standard w/ no basis",
        "clearing w/ basis",
        "clearing w/ no basis",
        "compression",
    ]

    for flag, label in zip(flags, labels):

        text = []
        tcon = []
        tred = []

        X = freudenthal_grid_light(m, n)
        for img in tqdm(data):
            irev = img

            t0 = time.monotonic()
            vals, imap = bats.lower_star_filtration(X, irev.flatten())
            t1 = time.monotonic()
            text.append(t1 - t0)

            t0 = time.monotonic()
            F = bats.FilteredLightSimplicialComplex(X, vals)
            t1 = time.monotonic()
            tcon.append(t1 - t0)

            t0 = time.monotonic()
            R = bats.reduce(F, bats.F2(), *flag)
            t1 = time.monotonic()
            tred.append(t1 - t0)

        print("\n\nBATS freudenthal {}".format(label))
        print("extension: {} sec.".format(np.mean(text)))
        print("construction: {} sec.".format(np.mean(tcon)))
        print("tred: {} sec.".format(np.mean(tred)))
        print("avg total: {} sec.".format(np.mean(text) + np.mean(tcon) + np.mean(tred)))


        text = []
        tcon = []
        tred = []

        X = cubical_grid(m,n)
        for img in tqdm(data):
            irev = img

            t0 = time.monotonic()
            vals = bats.lower_star_filtration(X, irev)
            t1 = time.monotonic()
            text.append(t1 - t0)

            t0 = time.monotonic()
            F = bats.FilteredCubicalComplex(X, vals)
            t1 = time.monotonic()
            tcon.append(t1 - t0)

            t0 = time.monotonic()
            R = bats.reduce(F, bats.F2(), *flag)
            t1 = time.monotonic()
            tred.append(t1 - t0)

        print("\n\nBATS cubical {}".format(label))
        print("extension: {} sec.".format(np.mean(text)))
        print("construction: {} sec.".format(np.mean(tcon)))
        print("tred: {} sec.".format(np.mean(tred)))
        print("avg total: {} sec.".format(np.mean(text) + np.mean(tcon) + np.mean(tred)))



    def bats_update_time(aimg, data, *flags):

        X = freudenthal_grid_light(m, n)
        text = []
        tupd = []

        # initialize on average image
        irev = aimg
        vals, imap = bats.lower_star_filtration(X, irev.flatten())
        F = bats.FilteredLightSimplicialComplex(X, vals)
        R = bats.reduce(F, bats.F2(), *flags)

        for img in tqdm(data):
            irev = img
            R = bats.reduce(F, bats.F2()) # reduce from scratch

            t0 = time.monotonic()
            vals, imap = bats.lower_star_filtration(X, irev.flatten())
            t1 = time.monotonic()
            text.append(t1 - t0)

            t0 = time.monotonic()
            R.update_filtration(vals)
            t1 = time.monotonic()
            tupd.append(t1 - t0)

        print("extension: {} sec.".format(np.mean(text)))
        print("update: {} sec.".format(np.mean(tupd)))
        print("avg total: {} sec.".format(np.mean(text) + np.mean(tupd)))



    print("\n\nBATS update average image:")
    aimg = np.mean(data, axis=0)
    bats_update_time(aimg, data)
    print("\nwith clearing basis")
    bats_update_time(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())

    # print("\n\nBATS update zero image:")
    # aimg = np.zeros((n,n))
    # bats_update_time(aimg, data)

    print("\n\nBATS cubical average image:")
    def time_BATS_cube(aimg, data, *flags):

        text = []
        tupd = []

        m, n = aimg.shape
        X = cubical_grid(m,n)

        # initialize on average image
        irev = aimg
        vals = bats.lower_star_filtration(X, irev)
        F = bats.FilteredCubicalComplex(X, vals)
        R = bats.reduce(F, bats.F2(), *flags)

        for img in tqdm(data):
            irev = img
            R = bats.reduce(F, bats.F2()) # reduce from scratch

            t0 = time.monotonic()
            vals = bats.lower_star_filtration(X, irev)
            t1 = time.monotonic()
            text.append(t1 - t0)

            t0 = time.monotonic()
            R.update_filtration(vals)
            t1 = time.monotonic()
            tupd.append(t1 - t0)

        print("extension: {} sec.".format(np.mean(text)))
        print("update: {} sec.".format(np.mean(tupd)))
        print("avg total: {} sec.".format(np.mean(text) + np.mean(tupd)))

    time_BATS_cube(aimg, data)
    print("\nwith clearing basis")
    time_BATS_cube(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())

    print("\n\ndionysus")
    def time_dionysus(data):

        tdion = []
        for img in data:
            irev = img
            t0 = time.monotonic()
            f_lower_star = dion.fill_freudenthal(irev)
            p = dion.homology_persistence(f_lower_star)
            dgms = dion.init_diagrams(p, f_lower_star)
            t1 = time.monotonic()
            tdion.append(t1 - t0)

        return tdion

    tdion = time_dionysus(data)
    print("dion avg: {} sec.".format(np.mean(tdion)))

    print("\n\nGudhi")
    def construct_gudhi_simplex_tree(img):
        # create in BATS to dump
        X = freudenthal_grid_light(*img.shape)
        vals, imap = bats.lower_star_filtration(X, img.flatten())
        F = bats.FilteredLightSimplicialComplex(X, vals)

        t0 = time.monotonic()
        GT = gd.SimplexTree()
        for k in range(F.maxdim() + 1):
            for s, v in zip(X.get_simplices(k), F.vals(k)):
                GT.insert(s, v)
        t1 = time.monotonic()
        return GT, t1 - t0

    def time_gudhi(data):

        tgd = []

        for img in data:
            irev = img

            GT, tcon = construct_gudhi_simplex_tree(irev)

            t0 = time.monotonic()
            diag = GT.persistence()
            t1 = time.monotonic()
            tgd.append(tcon + (t1 - t0))

        return tgd

    tgd = time_gudhi(data)
    print("gudhi avg: {} sec.".format(np.mean(tgd)))

    print("\n\nGudhi cubical")
    def time_gudhi_cubical(data):

        tgd = []

        for img in data:
            irev = img
            t0 = time.monotonic()
            cc = gd.CubicalComplex(
                dimensions = list(img.shape),
                top_dimensional_cells = img.flatten()
            )

            diag = cc.persistence()
            t1 = time.monotonic()
            tgd.append(t1 - t0)

        return tgd

    tgd = time_gudhi_cubical(data)
    print("gudhi avg: {} sec.".format(np.mean(tgd)))

    print("\n\nRipser lower star")
    def time_ripser(data):

        trips = []

        for img in data:
            irev = 255.0 - img
            t0 = time.monotonic()
            dgm = ripser.lower_star_img(img)
            t1 = time.monotonic()
            trips.append(t1 - t0)

        return trips

    tripser = time_ripser(data)
    print("ripser avg: {} sec.".format(np.mean(tripser)))

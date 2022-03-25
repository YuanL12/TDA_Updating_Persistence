import bats
import gudhi as gd
import dionysus as dion
import ripser
import keras
import matplotlib.pyplot as plt
from freudenthal import *
import time
import numpy as np
from tqdm import tqdm
import ripser

from keras.datasets import mnist

from kendall_tau import normalized_kt

def run_all():

    times = dict()
    lstr = 'MNIST'
    times[lstr] = dict()
    times[lstr]['Freudenthal'] = dict()
    times[lstr]['Cubical'] = dict()


    (train_X, train_y), (test_X, test_y) = mnist.load_data()
    print('X_train: ' + str(train_X.shape))
    print('Y_train: ' + str(train_y.shape))
    print('X_test:  '  + str(test_X.shape))
    print('Y_test:  '  + str(test_y.shape))

    img = train_X[0]
    img = 255 - img

    m, n = img.shape
    X = bats.LightFreudenthal(m, n)

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

    ps = R.persistence_pairs(0, False) +  R.persistence_pairs(1, False)
    for p in ps:
        if p.length() > 0:
            print(p)

    N = 1000
    data = train_X[:N]

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
        "compression w/ no basis",
    ]
    short_labels = [
        "BATS",
        "BATS(nb)",
        "BATS(cl)",
        "BATS(clb)",
        "BATS(co)"
    ]

    # for flag, label, sl in zip(flags, labels, short_labels):
    flag = (bats.standard_reduction_flag(), bats.clearing_flag())
    label = 'clearing w/ no basis'
    sl = 'BATS(cl)'

    text = []
    tcon = []
    tred = []
    X = bats.LightFreudenthal(m, n)
    for img in tqdm(data):
        irev = 255 - img

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
    times[sl + 'f'] = np.mean(text) + np.mean(tcon) + np.mean(tred)
    times[lstr]['Freudenthal']['BATS(c)'] = np.mean(text) + np.mean(tcon) + np.mean(tred)

    text = []
    tcon = []
    tred = []

    X = cubical_grid(m,n)
    for img in tqdm(data):
        irev = 255-img

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
    times[sl + 'c'] = np.mean(text) + np.mean(tcon) + np.mean(tred)
    times[lstr]['Cubical']['BATS(c)'] = np.mean(text) + np.mean(tcon) + np.mean(tred)


    def bats_update_time(aimg, data, *flags):

        X = bats.LightFreudenthal(m, n)
        text = []
        tupd = []
        ktds = []

        # initialize on average image
        irev = 255-aimg
        vals, imap = bats.lower_star_filtration(X, irev.flatten())
        F = bats.FilteredLightSimplicialComplex(X, vals)
        R = bats.reduce(F, bats.F2(), *flags)

        for img in tqdm(data):
            irev = 255-img
            R = bats.reduce(F, bats.F2()) # reduce from scratch

            t0 = time.monotonic()
            vals2, imap = bats.lower_star_filtration(X, irev.flatten())
            t1 = time.monotonic()
            text.append(t1 - t0)

            t0 = time.monotonic()
            R.update_filtration(vals2)
            t1 = time.monotonic()
            tupd.append(t1 - t0)

            ktds.append(normalized_kt(vals, vals2))



        print("extension: {} sec.".format(np.mean(text)))
        print("update: {} sec.".format(np.mean(tupd)))
        print("avg total: {} sec.".format(np.mean(text) + np.mean(tupd)))
        print("avg kt: {}".format(np.mean(ktds)))
        return np.mean(text) + np.mean(tupd), np.mean(ktds)



    print("\n\nBATS update average image:")
    aimg = np.mean(data, axis=0)
    times["BATS(u) avg"], times['kt avg'] = bats_update_time(aimg, data)
    print("\nwith clearing basis:")
    times["BATS(u,cl) avg"], _ = bats_update_time(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())

    print("\n\nBATS update zero image:")
    aimg = np.zeros((28,28))
    times["BATS(u) 0"], times['kt 0']= bats_update_time(aimg, data)
    print("\nwith clearing basis:")
    times["BATS(u,cl) 0"], _ = bats_update_time(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())


    print("\n\nBATS update random image:")
    aimg = np.random.rand(28,28)*255
    times["BATS(u) rnd"], times['kt rnd'] = bats_update_time(aimg, data)
    print("\nwith clearing basis:")
    times["BATS(u,cl) rnd"], _ = bats_update_time(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())


    print("\n\nBATS update selected image:")
    aimg = aimg = data[5]
    times["BATS(u) img"], times['kt img'] = bats_update_time(aimg, data)
    print("\nwith clearing basis:")
    times["BATS(u,cl) img"], _ = bats_update_time(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())

    times[lstr]['Freudenthal']['BATS(u,s)'], times[lstr]['Freudenthal']['dK'] = bats_update_time(aimg, data)
    times[lstr]['Freudenthal']['BATS(u,c)'], _ = bats_update_time(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())




    def time_BATS_cube(aimg, data, *flags):

        text = []
        tupd = []
        ktds = []

        m, n = aimg.shape
        X = cubical_grid(m,n)

        # initialize on average image
        irev = 255-aimg
        vals = bats.lower_star_filtration(X, irev)
        F = bats.FilteredCubicalComplex(X, vals)
        R = bats.reduce(F, bats.F2(), *flags)

        for img in tqdm(data):
            irev = 255-img
            R = bats.reduce(F, bats.F2()) # reduce from scratch

            t0 = time.monotonic()
            vals2 = bats.lower_star_filtration(X, irev)
            t1 = time.monotonic()
            text.append(t1 - t0)

            t0 = time.monotonic()
            R.update_filtration(vals2)
            t1 = time.monotonic()
            tupd.append(t1 - t0)

            ktds.append(normalized_kt(vals, vals2))

        print("extension: {} sec.".format(np.mean(text)))
        print("update: {} sec.".format(np.mean(tupd)))
        print("avg total: {} sec.".format(np.mean(text) + np.mean(tupd)))
        print("avg kt: {}".format(np.mean(ktds)))
        return np.mean(text) + np.mean(tupd), np.mean(ktds)

    print("\n\nBATS cubical average image:")
    times["BATS(u)c avg"], times['ktc avg'] = time_BATS_cube(aimg, data)
    print("\nwith clearing basis:")
    times["BATS(u,cl)c avg"], _ = time_BATS_cube(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())

    print("\n\nBATS cubical average image:")
    aimg = np.mean(data, axis=0)
    times["BATS(u)c 0"], times['ktc 0'] = time_BATS_cube(aimg, data)
    print("\nwith clearing basis:")
    times["BATS(u,cl) 0"], _ = time_BATS_cube(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())

    print("\n\nBATS cubical random image:")
    aimg = np.random.rand(28,28)*255
    times["BATS(u)c rnd"], times['ktc rnd'] = time_BATS_cube(aimg, data)
    print("\nwith clearing basis:")
    times["BATS(u,cl) rnd"], _ = time_BATS_cube(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())

    print("\n\nBATS cubical selected image:")
    aimg = aimg = data[5]
    times["BATS(u)c img"], times['ktc img'] = time_BATS_cube(aimg, data)
    print("\nwith clearing basis:")
    times["BATS(u,cl) img"], _ = time_BATS_cube(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())

    times[lstr]['Cubical']['BATS(u,s)'], times[lstr]['Cubical']['dK'] = time_BATS_cube(aimg, data)
    print("\nwith clearing basis")
    times[lstr]['Cubical']['BATS(u,c)'], _ = time_BATS_cube(aimg, data, bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag())


    print("\n\ndionysus")
    def time_dionysus(data):

        tdion = []
        for img in tqdm(data):
            irev = 255.0 - img
            t0 = time.monotonic()
            f_lower_star = dion.fill_freudenthal(irev)
            p = dion.homology_persistence(f_lower_star)
            dgms = dion.init_diagrams(p, f_lower_star)
            t1 = time.monotonic()
            tdion.append(t1 - t0)

        return tdion

    tdion = time_dionysus(data)
    print("dion avg: {} sec.".format(np.mean(tdion)))
    times['dion'] = np.mean(tdion)
    times[lstr]['Freudenthal']['Dionysus'] = np.mean(tdion)

    print("\n\nGudhi")
    def construct_gudhi_simplex_tree(img):
        # create in BATS to dump
        X = bats.LightFreudenthal(*img.shape)
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

        for img in tqdm(data):
            irev = 255.0 - img

            GT, tcon = construct_gudhi_simplex_tree(irev)

            t0 = time.monotonic()
            diag = GT.persistence()
            t1 = time.monotonic()
            tgd.append(tcon + (t1 - t0))

        return tgd

    tgd = time_gudhi(data)
    print("gudhi avg: {} sec.".format(np.mean(tgd)))
    times['gudhi(f)'] = np.mean(tgd)
    times[lstr]['Freudenthal']['Gudhi']  = np.mean(tgd)

    print("\n\nGudhi cubical")
    def time_gudhi_cubical(data):

        tgd = []

        for img in tqdm(data):
            irev = 255.0 - img
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
    times['gudhi(c)'] = np.mean(tgd)
    times[lstr]['Cubical']['Gudhi']  = np.mean(tgd)


    print("\n\nRipser lower star")
    def time_ripser(data):

        trips = []

        for img in tqdm(data):
            irev = 255.0 - img
            t0 = time.monotonic()
            dgm = ripser.lower_star_img(img)
            t1 = time.monotonic()
            trips.append(t1 - t0)

        return trips

    tripser = time_ripser(data)
    print("ripser avg: {} sec.".format(np.mean(tripser)))
    times['ripser'] = np.mean(tripser)
    times[lstr]['Freudenthal']['Ripser']  = np.mean(tripser)

    return times


if __name__ == "__main__":
    times = run_all()
    print(times)

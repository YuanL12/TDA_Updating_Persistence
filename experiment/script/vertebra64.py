import bats
import gudhi as gd
import dionysus as dion
# from freudenthal import levelset_cube
from freudenthal import levelset_cube_light, cubical_cube
import numpy as np
import time

from kendall_tau import normalized_kt

def run_all():

    times = dict()
    lstr = 'Vert-64'
    times[lstr] = dict()
    times[lstr]['Cubical'] = dict()


    # dimension of vertebra data
    dims = (512, 512, 512)

    # load vertebra data
    # A = np.fromfile('/project/bradnelson/update_persistence/levelset/data/vertebra8.raw', dtype=np.int8)
    A = np.fromfile('/home/brad/code/TDA_Updating_Persistence/ipynb/data/vertebra8.raw', dtype=np.int8)
    A = A.reshape(dims)
    img = np.array(A[::8,::8,::8], copy=True) # (each dimension is 512/8 = 64)

    # update from image plus noise
    img2 = img + 0.01 * np.random.randn(*img.shape)

    print("image shape: {}".format(img.shape))


    print("\n\nGudhi cubical")
    def time_gudhi_cubical(img):
        irev = img
        t0 = time.monotonic()
        cc = gd.CubicalComplex(
            dimensions = list(img.shape),
            top_dimensional_cells = img.flatten()
        )
        diag = cc.persistence()
        t1 = time.monotonic()

        return t1 - t0

    tgd = time_gudhi_cubical(img2)
    print("gudhi avg: {} sec.".format(np.mean(tgd)))
    times['gudhi(c)'] = np.mean(tgd)
    times[lstr]['Cubical']['Gudhi'] = np.mean(tgd)

    print("\n\nBATS cubical")
    dims = img.shape
    t0 = time.monotonic()
    X = bats.Cube(*dims)
    vals = bats.lower_star_filtration(X, img2)
    F = bats.FilteredCubicalComplex(X, vals)
    R = bats.reduce(F, bats.F2(), bats.standard_reduction_flag(), bats.clearing_flag())
    t1 = time.monotonic()
    times[lstr]['Cubical']['BATS(c)'] = t1 - t0
    print("BATS clearing: {} sec.".format(t1 - t0))


    def time_BATS_update_light_3D(img, img2, *flags):

        dims = img.shape
        X = bats.Cube(*dims)

        t0 = time.monotonic()
        vals = bats.lower_star_filtration(X, img)
        t1 = time.monotonic()
        print("time to extend: {} sec.".format(t1 - t0))

        t0 = time.monotonic()
        F = bats.FilteredCubicalComplex(X, vals)
        t1 = time.monotonic()
        print("time to construct: {} sec.".format(t1 - t0))

        t0 = time.monotonic()
        R = bats.reduce(F, bats.F2(), *flags)
        t1 = time.monotonic()
        print("time to reduce: {} sec.".format(t1 - t0))

        t0 = time.monotonic()
        vals2 = bats.lower_star_filtration(X, img2)
        R.update_filtration(vals2)
        t1 = time.monotonic()
        print("time to update: {} sec.".format(t1 - t0))
        tup = t1 - t0

        ktd = normalized_kt(vals, vals2)

        # t0 = time.monotonic()
        # vals = bats.lower_star_filtration(X, img2)
        # F = bats.FilteredCubicalComplex(X, vals)
        # R = bats.reduce(F, bats.F2(), *flags)
        # t1 = time.monotonic()
        # print("img2 from scratch: {} sec.".format(t1 - t0))
        # tstd = t1 - t0

        return tup, ktd

    print("\nBATS updates:")
    tup, ktd = time_BATS_update_light_3D(
        img,
        img2
    )
    times[lstr]['Cubical']['BATS(u,s)'] = tup
    times[lstr]['Cubical']['dK'] = ktd
    print("\nwith clearing")
    tup, ktd = time_BATS_update_light_3D(
        img,
        img2,
        bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag()
    )
    times[lstr]['Cubical']['BATS(u,c)'] = tup
    # times[lstr]['Cubical']['BATS(c)'] = tstd

    return times

if __name__ == "__main__":
    times = run_all()
    print(times)

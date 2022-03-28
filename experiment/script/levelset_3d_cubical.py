import bats
import gudhi as gd
import dionysus as dion
# from freudenthal import levelset_cube
from freudenthal import *
import numpy as np
import time

print("hello from levelset_3d_cubical.py")

# dimension of vertebra data
dims = (512, 512, 512)


t0 = time.monotonic()
# load vertebra data
A = np.fromfile('/project/bradnelson/yuan/levelset/data/vertebra8.raw', dtype=np.int8)
t1 = time.monotonic()
print("time to load data: {}".format(t1 - t0))
A = A.reshape(dims)
img = A[::8,::8,::8] # (each dimension is 512/8 = 64)
# img = A[::4,::4,::4] # (each dimension is 512/4 = 128)

# update from image plus noise
img = 1.0*img
img2 = img + 0.01 * np.random.randn(*img.shape)



print("image shape: {}".format(img.shape))
dims2 = img.shape

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

flags = [
    (),
    (bats.standard_reduction_flag(),),
    (bats.standard_reduction_flag(), bats.clearing_flag()),
    (bats.standard_reduction_flag(), bats.compression_flag()),
]
labels = [
    "standard w/ basis",
    "standard w/ no basis",
    "standard w/ clearing",
    "standard w/ compression",
]

t0 = time.monotonic()
X = cubical_cube(*dims2)
t1 = time.monotonic()
print("time to construct cubical complex: {}".format(t1 - t0))

for flag, label in zip(flags, labels):


    text = []
    tcon = []
    tred = []


    t0 = time.monotonic()
    vals = bats.lower_star_filtration(X, img2)
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


print("\n\nBATS cubical update image:")
def time_BATS_cube(img, img2):

    text = []
    tupd = []


    # initialize on average image
    vals = bats.lower_star_filtration(X, img)
    F = bats.FilteredCubicalComplex(X, vals)
    R = bats.reduce(F, bats.F2())

    t0 = time.monotonic()
    vals = bats.lower_star_filtration(X, img2)
    t1 = time.monotonic()
    text.append(t1 - t0)

    t0 = time.monotonic()
    R.update_filtration(vals)
    t1 = time.monotonic()
    tupd.append(t1 - t0)

    print("extension: {} sec.".format(np.mean(text)))
    print("update: {} sec.".format(np.mean(tupd)))
    print("avg total: {} sec.".format(np.mean(text) + np.mean(tupd)))

time_BATS_cube(img, img2)

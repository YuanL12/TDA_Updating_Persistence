import bats
import gudhi as gd
import dionysus as dion
# from freudenthal import levelset_cube
from freudenthal import levelset_cube_light
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
	A = np.fromfile('/project/bradnelson/yuan/levelset/data/vertebra8.raw', dtype=np.int8)
	A = A.reshape(dims)
	img = A[::8,::8,::8] # (each dimension is 512/8 = 64)

	# update from image plus noise
	img2 = img + 0.01 * np.random.randn(*img.shape)

	print("image shape: {}".format(img.shape))



	print("\n\nGudhi")
	def construct_gudhi_simplex_tree(img):
	    # create in BATS to dump
	    X = levelset_cube_light(*img.shape)
	    vals, imap = bats.lower_star_filtration(X, img.flatten())
	    F = bats.FilteredLightSimplicialComplex(X, vals)

	    t0 = time.monotonic()
	    GT = gd.SimplexTree()
	    for k in range(F.maxdim() + 1):
	        for s, v in zip(X.get_simplices(k), F.vals(k)):
	            GT.insert(s, v)
	    t1 = time.monotonic()
	    return GT, t1 - t0

	def time_gudhi(img):

	    GT, tcon = construct_gudhi_simplex_tree(img)

	    t0 = time.monotonic()
	    diag = GT.persistence()
	    t1 = time.monotonic()

	    return t1 - t0 + tcon

	# tgd = time_gudhi(img2)
	# print("gudhi avg: {} sec.".format(np.mean(tgd)))
	# times['gudhi'] = np.mean(tgd)



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

	print("\n\ndionysus")
	def time_dionysus(img):

	    tdion = []
	    t0 = time.monotonic()
	    f_lower_star = dion.fill_freudenthal(img)
	    p = dion.homology_persistence(f_lower_star)
	    dgms = dion.init_diagrams(p, f_lower_star)
	    t1 = time.monotonic()
	    tdion.append(t1 - t0)

	    return t1 - t0

	# tdion = time_dionysus(img2)
	# print("dion avg: {} sec.".format(np.mean(tdion)))
	# times['dionysus'] = np.mean(tdion)

	def time_BATS_update_light_3D(img, img2, *flags):

	    dims = img.shape
	    X = levelset_cube_light(*dims)

	    t0 = time.monotonic()
	    vals, imap = bats.lower_star_filtration(X, img.flatten())
	    t1 = time.monotonic()
	    print("time to extend: {} sec.".format(t1 - t0))

	    t0 = time.monotonic()
	    F = bats.FilteredLightSimplicialComplex(X, vals)
	    t1 = time.monotonic()
	    print("time to construct: {} sec.".format(t1 - t0))

	    t0 = time.monotonic()
	    R = bats.reduce(F, bats.F2(), *flags)
	    t1 = time.monotonic()
	    print("time to reduce: {} sec.".format(t1 - t0))

	    t0 = time.monotonic()
	    vals2, imap = bats.lower_star_filtration(X, img2.flatten())
	    R.update_filtration(vals2)
	    t1 = time.monotonic()
	    print("time to update: {} sec.".format(t1 - t0))
		tup = t1 - t0

		ktd = normalized_kt(vals, vals2)

	    t0 = time.monotonic()
	    vals, imap = bats.lower_star_filtration(X, img2.flatten())
	    F = bats.FilteredLightSimplicialComplex(X, vals)
	    R = bats.reduce(F, bats.F2(), *flags)
	    t1 = time.monotonic()
	    print("img2 from scratch: {} sec.".format(t1 - t0))
		tstd = t1 - t0

		return tsd, tup, ktd



	print("\nBATS:")
	tstd, tup, ktd = time_BATS_update_light_3D(
	    img,
	    img2
	)
	times[lstr]['Cubical']['BATS(u,s)'] = tup
	times[lstr]['Cubical']['dK'] = ktd
	print("\nwith clearing")
	tstd, tup, ktd = time_BATS_update_light_3D(
	    img,
	    img2,
		bats.standard_reduction_flag(), bats.clearing_flag(), bats.compute_basis_flag()
	)
	times[lstr]['Cubical']['BATS(u,c)'] = tup
	times[lstr]['Cubical']['BATS(c)'] = tstd

	return times

if __name__ == "__main__":
    times = run_all()
	print(times)

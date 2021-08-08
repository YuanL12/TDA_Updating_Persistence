### Reduction Algorithm
```C++
// perform reduction algorithm on a column matrix in-place
// apply change of basis to U
// invariant M * U = R
p2c_type reduce_matrix_standard(ColumnMatrix<TVec> &M, ColumnMatrix<TVec> &U) {
	std::cout << " go into the reduce_matrix_standard function in reduction.hpp" << std::endl;
	if (M.ncol() != U.ncol()) {throw std::runtime_error("Number of columns are not the same!");}

	// p2c_type pivot_to_col;
	p2c_type pivot_to_col(M.nrow(), bats::NO_IND);
	// create a temporary vector for use in axpys
	typename TVec::tmp_type tmp;

	// loop over columns
	for (size_t j = 0; j < M.ncol(); j++) {
		while(M[j].nnz() > 0) {
			// std::cout << j << " : ";
			// M[j].print_row();
			// piv is index-value nzpair
			auto piv = M[j].lastnz();
			// if (pivot_to_col.count(piv.ind) > 0) {
			if (pivot_to_col[piv.ind] != bats::NO_IND) { //Looking at if the pivot of column j has already been appeared in (pivot_to_col) previous column  is equivalent to look at the column before j with the same pivot as j
				// eliminate pivot
				size_t k = pivot_to_col[piv.ind];
				auto a = piv.val / M[k].lastnz().val;
				M[j].axpy(-a, M[k], tmp);
				U[j].axpy(-a, U[k], tmp); // update change of basis
			} else {
				// new pivot
				pivot_to_col[piv.ind] = j;
				break;
			}
		}
	}
	return pivot_to_col;
}
```

One thing worth noticing is the criterion for determining if there are columns before has the same pivot as the current column in the loop. 
BATs smartly uses `pivot_to_col` to store where each row appear as a pivot of a column, and so, for instance, `pivot_to_col[1]` will be the index of the column in which 1 (the index of the 2nd row) is the column's pivot. And, the index of a row does not appear as a pivot of any column, `pivot_to_col` will store the value as the maximum of `size_t`.
During the reduction loop, 
- if the `pivot_to_col[piv.ind] == bats::NO_IND`, it means the pivot of the current column is still unused by previous column and so can be stored by `pivot_to_col[piv.ind] = j;`.
- if the `pivot_to_col[piv.ind] != bats::NO_IND`, it means the pivot of the current column has already been used by some previous column and so we should do the reduction! Also, `pivot_to_col[piv.ind]` will tell which column is using the pivot.

##### Key insight of Reduction algorithm
The way of understanding the Reduction algorithm is interesting and provides several applications.
First of all, it is performing reduction algorithm on a column matrix in-place and keeping the invariant $M * U = R$: we initializing M as R and U as an identity matrix, and then performing reduction on M and recording the reduciton operations in U, until M becomes reduced.

But, in fact, since U is playing the role of recording column operations, U are not needed to assumed or initialized as an identity matrix. Thus, the only thing we only need to care is the invariant: $R = M * U$. One application is: 
If there is another upper triangle matrix U_1 and we now want to reduce $M' = M * U_1^{-1}$, we could just pass $U_1$ as the initializer of $U$ in the reduction algorithm. And in the end of the reduction, $RU^{-1} = M'$


#### Memebers of R
ReducedChainComplex (RCC) has four members:
```C++
std::vector<MT> upper_triangle_mat = RCC.U; // basis matrices
std::vector<MT> reduced_mat = RCC.R; // reduced matrices
std::vector<std::vector<size_t>> homology_basis_indices = RCC.I;
std::vector<p2c_type> pivot_to_column= RCC.p2c;
```

- The first two (upper triangle matrices U and reduced matrices R) are just factorization results from reduction algorithm. 
- `I` is the indices of homology revealing basis, which is obtained by: for each dimension k, finding every index of the zero column in `R[k]`, which is also not a pivot for the columns in `R[k+1]`.
- `p2c` is the "pivot to column" map: for each dimension k, its elements stores the index of the column, if row is a pivot of the column, otherwise will the maximum value of `size_t`.
  For instance, `p2c[1][2]` will be, in `R[1]`, the index of the column, in which 2 (the index of the 3rd row) is the column's pivot.


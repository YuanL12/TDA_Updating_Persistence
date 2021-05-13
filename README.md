# TDA_Updating_Persistence
This is a project about updating_persistence mentored by Dr. Nelson

## Updating
Since BATS is a git submodule, to you may need to use
```
git submodule update --remote
```
and
```
git pull --recurse-submodules
```

## Demo file
go to demo
```
make
```

this should create a file called "hello.out".  Try running it.
```
./hello.out
```


## Getting Started
### Matrix
If you are new to BATs, please go to `matrix_tutorial.cpp` to see the implementation of matrices.
### Persistence Homology
Then go to `persistence_tutorial.cpp` to see the implementation of persistence homology.
#### Example
The format of persistence pair:
\<dimension> : (\<birth_filtration_value>, \<death_filtration_value>) \<birth_index,death_index>
, where birth_index and death_index are generally the index of a simplex at dimension k(creates the homology) and the index of a simplex at dimension k+1(destroies the homology). 
```
persistence pair at dim 0
0 : (2,inf) <0,-1>
0 : (3,3) <1,0>
0 : (5,5) <2,1>

persistence pair at dim 1
1 : (5,6) <2,0>
```

### Permutation in BATs

#### Notaion
Notice that there are several ways to represent a permutation.
1. **Two-line notation**: one-to-one mapping from the first line to the second line.
```math
\sigma=\left(\begin{array}{lllll}
1 & 2 & 3 & 4 & 5 \\
2 & 5 & 4 & 3 & 1
\end{array}\right)
```
2. **One line notation**: only using the second line in a Two-line notation.
```math
(2\ 5\ 4\ 3\ 1)
```
3. **Cycle notation**: normal notation used in abstract algebra
```math
\left(\begin{array}{lllll}
1 & 2 & 3 & 4 & 5 \\
2 & 5 & 4 & 3 & 1
\end{array}\right)=(125)(34)=(34)(125)=(34)(512) .
```
4. **List of new indices**: storing the new indices of each element
```math
\{ 5\ 1\ 4\ 3\ 2\}
```
,which means, after permutation,
a) the index of the first element of a vector should be 5
b) the index of the second element of a vector should be 1
...

#### Connection
In BATs, a permutation of a vector takes the **One line notation**. To understand the implementation of the permutation function, 
```C++
// apply inverse permutation in-place
// O(nz log nz) where nz is #non-zeros
void ipermute(const std::vector<size_t>  &perm) {
for (size_t i = 0; i < nnz(); i++) {
indval[i].ind = perm[indval[i].ind];
}
sort();
}
```

there is a notable connection between **One line notation** and **List of new indices**:

Let $v$ is a list of new indices, then the permutation induced by $v$ is the inverse of $\sigma(v)$, where $\sigma(v)$ is the permutation with the same order of elements in $v$.

For example, the permutation induced by $\{ 5\ 1\ 4\ 3\ 2\}$ is equal to $(5\ 1\ 4\ 3\ 2)^{-1}$. 

## Updating process
Orginal process for computing PH
```C++
 using FT = ModP<int, 2>; // Field type Mod 2
 std::vector<double> f0 = {0.0, 0.1, 0.2}; // Filtration value on each vertex
 vals = lower_star_filtration(X, f0); // Filtration value on each simplex
 auto F = bats::Filtration(X, vals); // Build Filtration
 auto C = bats::Chain(F, FT()); //Build FilteredChainComplex
 auto R = bats::Reduce(C); //Build ReducedFilteredChainComplex
```

There are several options to update and we will first list several structs/classes used in the updating process.
Before we start, in all options, we assume that the filtration has already been established by:
`auto F = bats::Filtration(X, vals);`, where X is a simplicial complex and vals is its filtration value of each simplex. 

1. 创造一个Filtered Chain Complex：
   `auto C = bats::Chain(F, FT());`

   It will return a FilteredChainComplex, which has a member ChainComplex.
   Notice that, when initialling C, the constructor of FilteredChainComplex will
   1.1 自动根据filtration值排序,并存储到一个member perm
   1.2 存储`ChainComplex<MT> C = F.complex()`的Chain Complex作为一个member

2. 创造一个Reduced Filtered Chain Complex，即，开始Reduction Algorithm：
   `auto R = bats::Reduce(C);`

   It will return a ReducedFilteredChainComplex, which has a member ReducedChainComplex.
   Notice that, when initializing R, the constructor of ReducedFilteredChainComplex will 
   1.1 存储FilteredChainComplex的member：permutation, filtration value and chain complex
   1.2 初始化一个`ReducedChainComplex<MT> RC = C.complex()`的ReducedChainComplex

3. 通过ReducedChainComplex进行Reduction Algorithm:
    ```C++
   	// compute reduced chain complex from chain complex
	ReducedChainComplex(const ChainComplex<MT> &C) {
		size_t dmax = C.maxdim() + 1;
		//dim = C.dim;
		U.resize(dmax);
		R.resize(dmax);
		I.resize(dmax); //  a set of indices that can be used for a homology revealing basis.
		p2c.resize(dmax); // store the pivot of each cloumn 

		// TODO: can parallelize this
		for (size_t k = 0; k < dmax; k++) {
			U[k] = MT::identity(C.dim(k));
			R[k] = C.boundary[k];
			p2c[k] = reduce_matrix(R[k], U[k]);
		}

		set_indices();

	}
    ```

### Updating persistence comparison
Option 1. Updating filtration value in FilteredChainComplex 
i.e., updating C directly by 
```C++
 C.update_filtration(vals);
 R = bats::Reduce(C);
 R.print_summary()
```

Option 2. Updating filtration value in ReducedFilteredChainComplex 
i.e., updating R by
```C++
 R.update_filtration(vals);
 R.print_summary();
```
Next, in the struct function `update_filtration` of ReducedFilteredChainComplex in filtered_basis.hpp, `ReducedChainComplex<MT> RC;` is updated by `RC.permute_basis(iperm);`.

Now, go to the `permute_basis` function in basis.hpp. There are two main steps: one is permute rows and column U[k] and R[k] 


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

```math
a^2+b^2=c^2
```



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



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
$$
\sigma=\left(\begin{array}{lllll}
1 & 2 & 3 & 4 & 5 \\
2 & 5 & 4 & 3 & 1
\end{array}\right)
$$
2. **One line notation**: only using the second line in a Two-line notation.
$$(2\ 5\ 4\ 3\ 1)$$
3. **Cycle notation**: normal noation used in abstract algebra
$$
\left(\begin{array}{lllll}
1 & 2 & 3 & 4 & 5 \\
2 & 5 & 4 & 3 & 1
\end{array}\right)=(125)(34)=(34)(125)=(34)(512) .
$$
4. **List of new indices**: storing the new indices of each element
$$\{ 5\ 1\ 4\ 3\ 2\} $$
,which means, after permutation,
a) the index of the first element of a vector should be 5
b) the index of the second element of a vector should be 1
...

#### Connection
In BATs, a permutation of a vector takes the **One line notation**. To understand the implemenation of the permutation function, 
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
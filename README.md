# TDA_Updating_Persistence
This is a project about updating_persistence mentored by Dr. Nelson

## Setup
If you want to run C++ files, please use for submodule
```
git clone --recursive git@github.com:YuanL12/TDA_Updating_Persistence.git
```

If you only want to run python files, please follow the instruction in https://bnels.github.io/BATS.py/#/installation.


## Get Started 
### Python
The Python API in provided in BATS.py, which should be the easiest way to implement our functions.
We recommend try on our .ipynb files in `ipynb/Rips Radius.ipynb` for Rips filtration and `ipynb/Optimizing Levelset.ipynb` for Levelset filtration.

### C++
#### Setup
Since BATS is a submodule, you might need to 
In order to see how to compute PH, go to `demo/cpp/persistence_demo.cpp`, and then 
```Terminal
make persistence_demo.out
```

, this should create a file called "persistence_demo.out". Try running it.
``` Terminal
./persistence_demo.out
```

If you want to update on Rips Filtration then go to `demo/cpp/update_rips.cpp` and lower star Filtration in `demo/cpp/update_lower_star.cpp`.

### datasets
We appreciate the datasets provided on the internet, but due the file sizes, we are unable to upload them fully on Github. We list the sources of them below, if you want to try on real-world datasets:

1. The Stanford 3D Scanning Rep ository. http://graphics.stanford.edu/data/3Dscanrep
2. Volvis rep ository (archived). https://web.archive.org/web/20150307144939/http://volvis.org/
3. Datasets used in `A Roadmap for the Computation of Persistent Homology': https://github.com/n-otter/PH-roadmap/tree/master/data_sets (Otter, N., Porter, M. A., Tillmann, U., Grindrod, P., and Harrington, H. A. A roadmap for the computation of p ersistent homology. EPJ Data Science 6, 1 (2017))

## BATS(C++) Introduction
BATS is a git submodule, which includes C++ implementations of computational topology. If you want to develop it or have a deep understand of it, then read this section. For an algorithm analysis of reduction algorithm, you can see `demo/cpp/intro.md`.

### Matrix (Optional)
In order to create a matrix, go to `demo/cpp/matrix_demo.cpp` to see the implementation of matrices.


### Persistence homology
Standard process for computing PH
```C++ 
 auto F = bats::Filtration(X, vals); // Build a Filtration
 auto C = bats::Chain(F, FT()); //Build a Filtered Chain Complex
 auto R = bats::Reduce(C); //Build a Reduced Filtered Chain Complex
```

Now you are able to see the persistence pairs by:

```C++
std::cout << "\npersistence pair at dim 0" << std::endl;
for (auto& p: R.persistence_pairs(0)) 
    {std::cout << p.str() << std::endl;
}

std::cout << "\npersistence pair at dim 1" << std::endl;
for (auto& p: R.persistence_pairs(1)) {
    std::cout << p.str() << std::endl;
}
```
, and it will return 
```
persistence pair at dim 0
0 : (2,inf) <0,-1>
0 : (3,3) <1,0>
0 : (5,5) <2,1>

persistence pair at dim 1
1 : (5,6) <2,0>
```
Each line is an n-dimensional persistence pair in the format of 
```
<dimension> : (<birth_filtration_value>, <death_filtration_value>) <birth_index_of_n_simplex, death_index_n_plus_1_simplex>
```
, where birth_index and death_index are generally the index of a simplex at dimension k(creates the homology) and the index of a simplex at dimension k+1(destroies the homology). 

### Updating persistence
There are several options to update.
#### Fixed-sized filtration 
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

#### General filtration
If the complex of a filtration has been sorted by its filration values (which is the normal case), then, as shown in `demo/cpp/update_rips.cpp`,

```C++
// update Filtered Chain Complex
auto UI = bats::Update_info(F_X, F_Y); // build information
FCC.update_filtration_general(UI);
```
or 
```C++
// update Reduced Filtered Chain Complex
auto UI = bats::Update_info(F_X, F_Y); // build information
RFCC.update_filtration_general(UI);
```

For a general filtration whose complex is built in an arbitrary order, then you can also check `demo/cpp/update_general.cpp`.


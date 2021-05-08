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

```
persistence pair at dim 0
0 : (2,inf) <0,-1>
0 : (3,3) <1,0>
0 : (5,5) <2,1>

persistence pair at dim 1
1 : (5,6) <2,0>
```

#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>
#include <plu.hpp>

#define FT ModP<int, 2>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;

int main (int argc, char* argv[]) {
   
	size_t d = 2; // dimension of Euclidean Space
	size_t n = 10; // number of points to sample

	// maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 0.4);

    auto X = bats::sample_sphere<double>(d, n);
    // auto X_data = X.data;
    // X_data.print();
	auto dist = bats::Euclidean(); // metric

    /*
    Now let's exend a filtration, 
    */
    auto F = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
    print_filtration_info(F);

    /*
    Now let's reduce
    */
    auto FCC = bats::Chain(F, FT()); //FilteredChainComplex
    // auto RFCC = bats::Reduce(FCC); //ReducedFilteredChainComplex
    // std::cout << "Reduced matrix summary is" << std::endl;
    // RFCC.print_summary();
}
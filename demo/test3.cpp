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
    
    std::vector<std::vector<std::vector<size_t>>> Splx_1; // simplicial complex
    std::vector<std::vector<std::vector<size_t>>> Splx_2;
    // std::vector<std::vector<std::vector<size_t>>> Splx_3;
	size_t d = 2; // dimension of Euclidean Space
	size_t n = 10; // number of points to sample

	// maximum simplex dimension
    size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 0.4);

    auto X = bats::sample_sphere<double>(d, n);
    // X.data.print();

	auto dist = bats::Euclidean(); // metric


    auto F_X = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
    print_filtration_info(F_X);
    for(size_t i = 0; i <= F_X.maxdim(); i++) {
        auto A = F_X.complex().get_simplices(i);
        Splx_1.emplace_back(A);
    }

    auto Y = bats::sample_sphere<double>(d, 13);
    auto F_Y = bats::RipsFiltration<CpxT>(Y, dist, rmax, maxdim);
    print_filtration_info(F_Y);
    for(size_t i = 0; i <= F_Y.maxdim(); i++) {
        auto A = F_Y.complex().get_simplices(i);
        Splx_2.emplace_back(A);
    }        
    

    // auto Z = bats::sample_sphere<double>(d, 11);
    // auto F_Z = bats::RipsFiltration<CpxT>(Z, dist, rmax, maxdim);
    // print_filtration_info(F_Z);
    // for(size_t i = 0; i <= F_Z.maxdim(); i++) {
    //     auto A = F_Z.complex().get_simplices(i);
    //     Splx_3.emplace_back(A);
    // }  


    auto UI = Updating_filtration_info(Splx_1, Splx_2);
    for(auto info : UI.information){
        info.print_summary();
    }

    // compare complex 2 and complex 3
    // UI = Updating_filtration_info(Splx_2, Splx_3);
    // for(auto info : UI.information){
    //     info.print_summary();
    // }


    /*
    We first compute PH by original order
    */
    auto FCC = bats::Chain(F_X, FT()); //FilteredChainComplex
    auto RFCC = bats::Reduce(FCC); //ReducedFilteredChainComplex
    // std::cout << "Reduced matrix summary is" << std::endl;
    // RFCC.print_summary();

    /*
    Updating filtration option 1: rebuild FCC and RFCC
    */
    // {
    //     auto FCC2 = bats::Chain(F_Y, FT()); //FilteredChainComplex
    //     auto RFCC2 = bats::Reduce(FCC2); //ReducedFilteredChainComplex
    //     std::cout << "Reduced matrix summary is" << std::endl;
    //     RFCC2.print_summary();
    // }

    /*
    Updating filtration option 2: update FCC and rebuild RFCC
    */
    // {
    //     // FCC.update_filtration(UI, F_Y);
    //     auto RFCC2 = bats::Reduce(FCC); //ReducedFilteredChainComplex
    //     std::cout << "Reduced matrix summary is" << std::endl;
    //     RFCC2.print_summary();
    // }
}
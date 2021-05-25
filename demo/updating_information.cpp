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



int main () {
   

    std::vector<std::vector<std::vector<size_t>>> Splx_1;
    std::vector<std::vector<std::vector<size_t>>> Splx_2;

    std::vector<std::vector<size_t>> A = {{1},{0},{2},{4},{3}};
    std::vector<std::vector<size_t>> B = {{0},{1},{5},{3},{4},{2}};
    Splx_1.emplace_back(A);
    Splx_2.emplace_back(B);

    A = {{0, 2}, {3, 5}, {2, 3},{2, 4},{3,4}};
    B = {{2, 4}, {0,1}, {0, 2}, {1,2},{1,4},{3, 5},{0, 4}};
    Splx_1.emplace_back(A);
    Splx_2.emplace_back(B);

    A = {{2, 3, 4}};
    B = {{0, 1, 4}};
    Splx_1.emplace_back(A);
    Splx_2.emplace_back(B);


    auto UI = Updating_filtration_info(Splx_1, Splx_2);
    for(auto info : UI.information){
        info.print_summary();
    }

}
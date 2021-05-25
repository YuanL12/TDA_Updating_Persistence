#include <bats.hpp>
#include <util/io.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include "plu.hpp"

#define FT ModP<int, 2>
#define VT SparseVector<FT>
#define MT ColumnMatrix<VT>

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;

int main () {
    /*
    First we build a simplicial complex which can be used to extend filtrations
    */
    CpxT X(4,2);
    
    std::vector<size_t> s;
    
    s = {0}; X.add(s);
    s = {1}; X.add(s);
    s = {2}; X.add(s);
    s = {3}; X.add(s);

    s = {2,3}; X.add(s);
    s = {0,1}; X.add(s);
    s = {1,2}; X.add(s);
    s = {0,2}; X.add(s);
    s = {0,3}; X.add(s);

    s = {0,1,2}; X.add(s);
    X.print_summary();

    CpxT Y(3,2);

    s = {0}; Y.add(s);
    s = {1}; Y.add(s);
    s = {2}; Y.add(s);
    s = {0,1}; Y.add(s);
    s = {1,2}; Y.add(s);
    s = {0,2}; Y.add(s);
    s = {0,1,2}; Y.add(s);

    Y.print_summary();
    
    /*
    Now let's exend a filtration
    */
    std::vector<double> f0 = {0.0, 0.2, 0.1, 0.3};

    
    /*
    The above is wrapped by the lower_star_filtration function
    */
    auto vals = lower_star_filtration(X, f0); // filtration values 是根据 X 中 simplex 的排列顺序来设定的
    std::cout << "\nold filtration values are:" << std::endl;
    for (auto& valsk: vals) {
        for (auto& v: valsk) {
            std::cout << v << ", ";
        }
        std::cout << "\n";
    }
    
    /*
    Now let's build a filtration and reduce
    */
    auto F = bats::Filtration(X, vals);
    auto FCC = bats::Chain(F, FT());
    std::cout << "sorted boundary matrices:" << std::endl;
    for(auto & M: FCC.complex().boundary){
        M.print();
    }
    // auto R = bats::Reduce(C);
    // R.print_summary();
    //
    // for (auto& p: R.persistence_pairs(0)) {
    //     std::cout << p.str() << std::endl;
    // }
    //
    /*
    If we want to change the filtration we have a variety of options.
    First is to update the filtration on the FilteredChainComplex
    Then, we reduce the complex again
    */
    // update filtration
    f0 = {1.1, 1.0, 1.2};

    vals = lower_star_filtration(Y, f0);
    std::cout << "\nnew filtration value is" << std::endl;
    for (auto& valsk: vals) {
        for (auto& v: valsk) {
            std::cout << v << ", ";
        }
        std::cout << "\n";
    }

    auto F_2 = bats::Filtration(Y, vals);
    FCC = bats::Chain(F_2, FT());
    std::cout << "sorted boundary matrices:" << std::endl;
    for(auto & M: FCC.complex().boundary){
        M.print();
    }

    // std::cout << "\ndive into update_filtration function" << std::endl;
    // auto perms = bats::filtration_sortperm(vals);
    // print_2D_vectors(perms);

    // What we want? 找到两个permutation 直接的比较 (pass)
    // 方案1 从boundary matrix 来看差异，一列一列，
    // 方案2 从filtration value直接看差异

    // 已有：将filtration value 变为升序排列的 permutation
    /*
    第一步，先把两个filtration 都升序排列
    第二步，loop over each column of the original matrix, for the i-th col, 
    */


    // R = bats::Reduce(C);
    // R.print_summary();
    //
    // for (auto& p: R.persistence_pairs(0)) {
    //     std::cout << p.str() << std::endl;
    // }
    //
    // /*
    // The second option is to update the ReducedFilteredChainComplex
    // In some situations, this may be faster
    // */
    // // update filtration again
    // f0 = {0.0, 0.1, 0.2};
    // vals = lower_star_filtration(X, f0);
    // R.update_filtration(vals);
    // R.print_summary();
    //
    // for (auto& p: R.persistence_pairs(0)) {
    //     std::cout << p.str() << std::endl;
    // }

    /*
    The final option would be to just call Chain() again to build
    the updated chain complex.
    */
}
#include <bats.hpp>
#include <vector>
#include <iostream>
#include <plu.hpp>

using FT = ModP<int, 2>;

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;

int main() {


    /*
    First we build a simplicial complex which can be used to extend filtrations
    */
    CpxT X(3,2);
    
    std::vector<size_t> s;
    
    s = {0}; X.add(s);
    s = {1}; X.add(s);
    s = {2}; X.add(s);
    s = {0,1}; X.add(s);
    s = {1,2}; X.add(s);
    s = {0,2}; X.add(s);
    s = {0,1,2}; X.add(s);
    
    X.print_summary();
    
    /*
    Now let's exend a filtration
    */
    std::vector<double> f0 = {0.0, 0.1, 0.2};
    // lower star filtration
    auto vals = lower_star_filtration(X, f0);
    std::cout << "\nfiltration values are:" << std::endl;
    print_2D_vectors(vals);
    
    auto perms = bats::filtration_sortperm(vals);
    std::cout << "\nperms is:" << std::endl;
    print_2D_vectors(perms);


    /*
    Now let's build a filtration and reduce
    */
    auto F = bats::Filtration(X, vals); // Filtration
    auto C = bats::Chain(F, FT()); //FilteredChainComplex
    auto R = bats::Reduce(C); //ReducedFilteredChainComplex
    R.print_summary();
    
    for (auto& p: R.persistence_pairs(0)) {
        std::cout << p.str() << std::endl;
    }
    
    /*
    If we want to change the filtration we have a variety of options.
    First is to update the filtration on the FilteredChainComplex
    Then, we reduce the complex again
    */
    // update filtration
    f0 = {1.1, 1.0, 1.2};
    vals = lower_star_filtration(X, f0);
    C.update_filtration(vals);
    R = bats::Reduce(C);
    R.print_summary();
    
    for (auto& p: R.persistence_pairs(0)) {
        std::cout << p.str() << std::endl;
    }
    
    return EXIT_SUCCESS;
}

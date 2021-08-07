#include <bats.hpp>
#include <vector>
#include <iostream>
#include <plu.hpp>
#include <vector>
#include <iostream>
#include <math.h>
#include <chrono>
#include <random>

using FT = ModP<int, 2>;

using CpxT = bats::LightSimplicialComplex<size_t, std::unordered_map<size_t, size_t>>;
// using CpxT = bats::SimplicialComplex;

int main() {


    /*
    First we build a simplicial complex which can be used to extend filtrations
    */
    CpxT X(7,2);
    
    std::vector<size_t> s;
    
    s = {0,1,2}; X.add_recursive(s);
    s = {0,1,3}; X.add_recursive(s);
    s = {0,1,4}; X.add_recursive(s);
    s = {0,1,5}; X.add_recursive(s);
    s = {0,1,6}; X.add_recursive(s);
    s = {0,2,3}; X.add_recursive(s);
    s = {0,2,4}; X.add_recursive(s);
    s = {0,2,5}; X.add_recursive(s);
    s = {0,2,6}; X.add_recursive(s);
    s = {0,3,4}; X.add_recursive(s);
    s = {0,3,5}; X.add_recursive(s);
    s = {0,3,6}; X.add_recursive(s);
    s = {0,4,5}; X.add_recursive(s);
    s = {0,4,6}; X.add_recursive(s);
    s = {0,5,6}; X.add_recursive(s);
    X.print_summary();
    
    /*
    Now let's exend a filtration
    */
    std::vector<double> f0 = {0.1, 0.2, 0.0, 0.5, 0.4, 0.8, 0.0};
    
    /*
    The above is wrapped by the lower_star_filtration function
    */
    auto vals = lower_star_filtration(X, f0);
    std::cout << "\nfiltration value is:" << std::endl;
    // for (auto& valsk: vals) {
    //     for (auto& v: valsk) {
    //         std::cout << v << ", ";
    //     }
    //     std::cout << "\n";
    // }
    
    /*
    Now let's build a filtration and reduce
    */
    auto F = bats::Filtration(X, vals); // Filtration
    auto FCC = bats::Chain(F, FT()); //FilteredChainComplex
    auto CC = FCC.C;
    auto boundary_mat = CC.boundary;

    auto RFCC = bats::Reduce(FCC); //ReducedFilteredChainComplex
    RFCC.print_summary();

    auto RCC = RFCC.RC; // ReducedChainComplex has four members:
	auto upper_triangle_mat = RCC.U; // basis matrices
	auto reduced_mat_mat = RCC.R; // reduced matrices
    // I is the set of indices for zero columns which do not appear as pivots one dimension higher.  
    // This is gives a set of indices that can be used for a homology revealing basis.
	std::vector<std::vector<size_t>> homology_basis_indices = RCC.I;
	std::vector<p2c_type> pivot_to_column= RCC.p2c;

    {
        /*
        If we want to change the filtration we have a variety of options.
        The naive way is to rebuild the filtration
        */
        // update filtration
        f0 = {0.0, 0.1, 0.2, 0.2, 0.0, 0.5, 0.4};
        vals = lower_star_filtration(X, f0);
        auto F = bats::Filtration(X, vals); // Filtration
        FCC = bats::Chain(F, FT()); //FilteredChainComplex

        RFCC = bats::Reduce(FCC); //ReducedFilteredChainComplex
        RFCC.print_summary();


        std::cout << "\nupdate the filtration by rebuiding filtration" << std::endl;
        std::cout << "persistence_pairs is" << std::endl;
        for (auto& p: RFCC.persistence_pairs(0)) {
            std::cout << p.str() << std::endl;
        }
        std::cout << "Reduced matrices are" << std::endl;
        auto RCC = RFCC.RC; // ReducedChainComplex has four members
	    auto reduced_mat_mat = RCC.R; // reduced matrices
        for (const auto& mat: reduced_mat_mat){
            mat.print();
        }
        std::cout << "\n" << std::endl;
    }
    
    {
        /*
        If we want to change the filtration we have a variety of options.
        Second is to update the filtration on the FilteredChainComplex
        Then, we reduce the complex again
        */
        // update filtration
        f0 = {0.0, 0.1, 0.2, 0.2, 0.0, 0.5, 0.4};
        vals = lower_star_filtration(X, f0);
        FCC.update_filtration(vals);
        RFCC = bats::Reduce(FCC);
        RFCC.print_summary();
        
        std::cout << "\nupdate the filtration on the FilteredChainComplex" << std::endl;
        std::cout << "persistence_pairs is" << std::endl;
        for (auto& p: RFCC.persistence_pairs(0)) {
            std::cout << p.str() << std::endl;
        }
        std::cout << "Reduced matrices are" << std::endl;
        auto RCC = RFCC.RC; // ReducedChainComplex has four members
	    auto reduced_mat_mat = RCC.R; // reduced matrices
        for (const auto& mat: reduced_mat_mat){
            mat.print();
        }
        std::cout << "\n" << std::endl;
    }
    
    {   /*
        The thrid option is to update the ReducedFilteredChainComplex
        In some situations, this may be faster
        */
        
        f0 = {0.0, 0.1, 0.2, 0.2, 0.0, 0.5, 0.4};
        vals = lower_star_filtration(X, f0);
        RFCC.update_filtration(vals);
        RFCC.print_summary();
        
        std::cout << "\nupdate the ReducedFilteredChainComplex" << std::endl;
        std::cout << "persistence_pairs is" << std::endl;
        for (auto& p: RFCC.persistence_pairs(0)) {
            std::cout << p.str() << std::endl;
        }
        std::cout << "Reduced matrices are" << std::endl;
        auto RCC = RFCC.RC; // ReducedChainComplex has four members
	    auto reduced_mat_mat = RCC.R; // reduced matrices
        for (const auto& mat: reduced_mat_mat){
            mat.print();
        }
        std::cout << "\n" << std::endl;
    }

    {
        /*
        The fourth option is to update the ReducedFilteredChainComplex by UQL
        */

        f0 = {0.0, 0.1, 0.2, 0.2, 0.0, 0.5, 0.4};
        vals = lower_star_filtration(X, f0);
        RFCC.update_filtration_by_UQL(vals);
        RFCC.print_summary();
        
        std::cout << "\nupdate the ReducedFilteredChainComplex by UQL" << std::endl;
        std::cout << "persistence_pairs is" << std::endl;
        for (auto& p: RFCC.persistence_pairs(0)) {
            std::cout << p.str() << std::endl;
        }
        std::cout << "Reduced matrices are" << std::endl;
        auto RCC = RFCC.RC; // ReducedChainComplex has four members
	    auto reduced_mat_mat = RCC.R; // reduced matrices
        for (const auto& mat: reduced_mat_mat){
            mat.print();
        }
    }

    return EXIT_SUCCESS;
}
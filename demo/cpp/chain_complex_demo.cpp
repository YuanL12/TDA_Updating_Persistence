#include "../BATS/include/bats.hpp"
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
    CpxT X(3,2); // 3 is the size of vertex set, 2 is the size of maximum simplex dimension
    std::vector<size_t> s;

    s = {0}; X.add(s);
    s = {1}; X.add(s);
    s = {2}; X.add(s);
    s = {0,1}; X.add(s);
    s = {1,2}; X.add(s);
    s = {0,2}; X.add(s);
    s = {0,1,2}; X.add(s);

    X.print_summary(); //summary of simpilcial complex X

    /*
    Now let's exend a filtration
    */
    std::vector<double> f0 = {0.0, 0.1, 0.2}; // filtration value on vertex
    // lower star filtration
    std::function<double(const std::vector<size_t>&)> filtfn = [&f0](
        const std::vector<size_t>& s
    ) -> double {
        return f0[*std::max_element(s.begin(), s.end(), [&f0](size_t i, size_t j) {return f0[i] < f0[j];})];
    };
    
    // filtration values are stored by the index of order of vertices not by their filtration values!!!
    //print_summary_of_filtration(X, filtfn);

    // std::cout << "\nextend_filtration" << std::endl;
    // store filtration vaules of all simplices into a vector of vectors
    auto vals = extend_filtration(X, filtfn);
    // for (auto& valsk: vals) { // each dimension
    //     for (auto& v: valsk) {
    //         std::cout << v << ", "; // each simplex's filtration value
    //     }
    //     std::cout << "\n";
    // }

    /*
    The above can be also achieved by one-step: using lower_star_filtration function
    */
    //std::cout << "\nold lower_star_filtration" << std::endl;
    vals = lower_star_filtration(X, f0);
    //print_2D_vectors(vals);

    /*
    Now let's build a filtration and reduce
    */
    auto F = bats::Filtration(X, vals);
    auto FCC = bats::Chain(F, FT()); // filtered chain complex has ChainComplex<MT> C as one of its member

    // boundary matrix of chain complex are also able to see
    std::cout << "\nAll boundary matices" << std::endl;
    for (auto boundary_mat : FCC.C.boundary) {
        boundary_mat.print();
    }

    // Todo: set a filtration value of a simplex manually

    auto R = bats::Reduce(FCC); // reduced matrix
    std::cout << "\nReduced matrix summary" << std::endl;
    R.print_summary();

    // if you want to see persistence pairs for H1/H0, use
    std::cout << "\nPH(0) is:" << std::endl;
    auto ps = R.persistence_pairs(0);
    for (auto p : ps) {
        if (p.death > p.birth)
            std::cout << p.str() << " " << p.death - p.birth << std::endl;
    }

    return EXIT_SUCCESS;
}
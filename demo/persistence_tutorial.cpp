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
    CpxT X(3,2); // 3 is the size of vertex set, 2 is the size of maximum simplex dimension

    std::vector<size_t> s; //store temporary simplices 

    s = {0}; X.add(s);
    s = {1}; X.add(s);
    s = {2}; X.add(s);
    s = {0,1}; X.add(s);
    s = {1,2}; X.add(s);
    s = {0,2}; X.add(s);
    s = {0,1,2}; X.add(s);

    //summary of simpilcial complex X
    std::cout << "\nsummary of simpilcial complex X" << std::endl;
    X.print_summary(); 

    /*
    Now let's exend a filtration
    */
    std::vector<double> f0 = {5.0, 2.0, 3.0}; // filtration value on each vertex
    
    // lower star filtration
    // store filtration vaules of all simplices into a vector of vectors
    auto vals = lower_star_filtration(X, f0);
    std::cout << "\nfiltration values on each simplices" << std::endl;
    print_2D_vectors(vals);

    /*
    Now let's build a filtration and reduce
    */
    auto F = bats::Filtration(X, vals);
    auto C = bats::Chain(F, FT());
    auto R = bats::Reduce(C);
    
    std::cout << "\noverall information of Reduced matrix" << std::endl;
    R.print_summary();

    std::cout << "\npersistence pair at dim 0" << std::endl;
    for (auto& p: R.persistence_pairs(0)) {
        std::cout << p.str() << std::endl;
    }

    return EXIT_SUCCESS;
}

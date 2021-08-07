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
    // filtration value on each simplex
    std::vector<std::vector<double>> vals = {{1,3,2},{5,6,4},{7}}; 
    /*
    Now let's build a filtration 
    */
    auto F = bats::Filtration(X, vals);
    // if  you want to change filtration value of a specific simplex 
    s = {0,1,3};
    F.add(3.0, s);
    print_filtration_info(F);

    std::cout << "\nsummary of simpilcial complex X" << std::endl;
    X.print_summary(); 
    return EXIT_SUCCESS;
}

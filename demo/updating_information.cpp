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


template<typename T>
void print_simplex (const T& perm) {
    std::cout << "{" ;
    for (auto v = perm.begin(); v < perm.end() - 1; v++ ) {
        std::cout << *v << ", ";
    }
    std::cout << *(perm.end()-1) << "} " ;
}

/*
Store the information of two lists of simplices, 
after set operations(intersection and complements),
which can be use to determine permutation/addition/deletion
*/
struct simplex_info
{
    std::vector<size_t> simplex;
    size_t index_in_A;
    size_t index_in_B;

    simplex_info(const std::vector<size_t> &splx, 
                        const size_t &ind_A, 
                        const size_t &ind_B): 
                        simplex(splx), 
                        index_in_A(ind_A),
                        index_in_B(ind_B){}
 
    void print(){
        print_simplex(simplex);
        std::cout << index_in_A << " "<< index_in_B << std::endl;
    }
};


/*
updating information at dimension k, passing into two lists of k-simplices
*/
struct Updating_filtration_info_k{
    std::vector<simplex_info> intersection_of_simplices;
    std::vector<simplex_info> deletion_of_simplices;
    std::vector<simplex_info> addition_of_simplices;
    size_t k; // dimension of the simplices

    Updating_filtration_info_k(const std::vector<std::vector<size_t>> &A, 
    const std::vector<std::vector<size_t>> &B){
        k = A[0].size() - 1;

        for (size_t i = 0; i < A.size(); i++){// find A \cap B and A - (A \cap B) 
            auto it = std::find(B.begin(), B.end(), A[i]);
            if (it != std::end(B)){ // permuation information, since A[i] is finded in B successfully
                auto index_it = std::distance(B.begin(), it); // find the index of finded simplex in B
                // std::cout << "B contains ";
                // print_simplex(*it);
                // std::cout << ", and its index in A is "<< i << ", and in B is " << index_it << std::endl;
                simplex_info inter_splx(*it, i, index_it);
                intersection_of_simplices.emplace_back(inter_splx);
            }
            else{ // deletion information, since A[i] is not finded in B
                // std::cout << "B does not contain ";
                // print_simplex(A[i]);
                // std::cout << "\n";
                simplex_info delete_splx(A[i], i, bats::NO_IND);
                deletion_of_simplices.emplace_back(delete_splx);
            } 
        }

        // next we find the addition information, i.e., find the elements in B but not in A \cap B
        for (size_t i = 0; i < B.size(); i++){
            bool find_or_not = false;
            for (auto iter = intersection_of_simplices.begin(); 
                iter < intersection_of_simplices.end(); iter++){
                    if (iter->simplex == B[i]) find_or_not = true;
            } //define 
            if (!find_or_not){
                simplex_info add_splx(B[i], bats::NO_IND, i);
                addition_of_simplices.emplace_back(add_splx);
            }
        }  

    }

    void print_summary(){
        std::cout << "\nsimplex lists of dimension "<< k << std::endl;
        std::cout << "print intersection" << std::endl;
        for (auto & v: intersection_of_simplices){
            v.print();
        }
        std::cout << "print deletion" << std::endl;
        for (auto & v: deletion_of_simplices){
            v.print();
        }
        std::cout << "print addition" << std::endl;
        for (auto & v: addition_of_simplices){
            v.print();
        }

    }
};

/*
Store the information after updating filtration values.
Need to pass into two soreted sets of simplicial complexes.
*/
template <typename simplicial_complex>
struct Updating_filtration_info
{
    std::vector<Updating_filtration_info_k> infomation;

    Updating_filtration_info(const simplicial_complex& splx1, 
    const simplicial_complex& splx2 ){
        for (size_t i=0; i<splx1.size(); i++){
            infomation.emplace_back(Updating_filtration_info_k(splx1[i], splx2[i]));
        }
    } 
};


template <class Filtration>
void print_filtration_info(const Filtration& F){

    auto F_complex = F.complex();
    std::cout << "\nRips Filtration Values are" << std::endl;
    auto filtration_vals = F.vals();

    for (size_t i = 0; i <= F.maxdim(); i++) {
        std::cout << F.ncells(i) << " cells in dim " << i << ":"<< std::endl;
        auto simplices_i = F_complex.get_simplices(i);
        for (size_t j = 0; j < F.ncells(i); j++){
            print_simplex(simplices_i[j]);
            std::cout << ": "<< filtration_vals[i][j] << std::endl;   
        }
    }
}

int main (int argc, char* argv[]) {
   
	// size_t d = 2; // dimension of Euclidean Space
	// size_t n = 10; // number of points to sample

	// // maximum simplex dimension
    // size_t maxdim = bats::util::io::parse_argv(argc, argv, "-maxdim", 3);
    // double rmax = bats::util::io::parse_argv(argc, argv, "-rmax", 0.4);

    // auto X = bats::sample_sphere<double>(d, n);
    // // auto m = X_data.m;
    // std::cout << "\ndata set is" << std::endl;
    // // X.data.print();
	// auto dist = bats::Euclidean(); // metric

    // /*
    // 如何构造complex的？
    // */
    // auto F = bats::RipsFiltration<CpxT>(X, dist, rmax, maxdim);
    // print_filtration_info(F);

    std::vector<std::vector<std::vector<size_t>>> Splx_1;
    std::vector<std::vector<std::vector<size_t>>> Splx_2;

    std::vector<std::vector<size_t>> A = {{0},{1},{2},{3},{4},{5},{6},{7}};
    std::vector<std::vector<size_t>> B = {{0},{1},{2},{3},{4},{5},{6}};
    Splx_1.emplace_back(A);
    Splx_2.emplace_back(B);

    A = {{1,6},{3,7},{1,5},{5,6}};
    B = {{2,3},{1,6},{5,6},{1,5},{2,4}};
    Splx_1.emplace_back(A);
    Splx_2.emplace_back(B);

    A = {{1,3,5},{2,3,6},{0,2,3}};
    B = {{1,2,4},{1,3,5}};
    Splx_1.emplace_back(A);
    Splx_2.emplace_back(B);

    auto UI = Updating_filtration_info(Splx_1, Splx_2);
    for(auto info : UI.infomation){
        info.print_summary();
    }

}
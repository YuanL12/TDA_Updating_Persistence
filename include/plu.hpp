/*
header file for PLU factorization
*/
#pragma once
#include <bats.hpp>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

template <typename CpxT, typename T>
void print_summary_of_filtration(const CpxT& X, std::function<T(const std::vector<size_t>&)>& filtfn){
    std::cout << "\nLet's see the filtration value on each dimension" << std::endl;
    for (size_t k = 0; k < X.maxdim() + 1; k++) {
        std::cout << "For dimension "<< k << std::endl;
        for (auto& s : X.get_simplices(k)) {
            std::cout << "simplex with index "<< X.find_idx(s) << ", filtration value f(s) is " << filtfn(s) << std::endl; 
        }
    }
}


//two ways to display 2D vectors
template<typename T>
void print_2D_vectors (const T& perms) {
//perms represents permutations
    for (auto& valsk: perms) {
        for (auto& v: valsk) {
            std::cout << v << " ";
        }
        std::cout << "\n";
    }
}
template<typename T>
void print_1D_vectors (const T& perm) {
    for (auto& v: perm) {
        std::cout << v << " ";
    }
    std::cout << "\n";
}

// Reduction algorithm of matrix permutation version 
// of updating persistence by passing old RU factorzation 
// and row and column permutations 
template <class TC>
void update_reduction(ColumnMatrix<TC> &R_k, ColumnMatrix<TC> &U_k, const ColumnMatrix<TC> & row_perm, const ColumnMatrix<TC> &column_perm){
	
	using MatT = ColumnMatrix<TC>;
	MatT P_k_minus_1 = row_perm;
	MatT P_k = column_perm;
	    
	auto F = UQL(P_k.T() * U_k);                                                
    MatT U = F.U;        
    MatT L = F.L;
    MatT Q = F.E;

    MatT R_k_prime = P_k_minus_1 * R_k * l_inv(L) * Q.T();
    MatT U_double_prime = MatT::identity(R_k_prime.ncol());
    bats::reduce_matrix(R_k_prime, U_double_prime);

    MatT U_k_prime = U * U_double_prime;

    R_k = R_k_prime;
    U_k = U_k_prime;    
}

// Reduction algorithm of list permutation version(maybe faster)
// for updating persistence by passing old RU factorzation 
// and row and column permutations (in the form of list permutation)
template <class TC>
void update_reduction(ColumnMatrix<TC> &R_k, ColumnMatrix<TC> &U_k, const std::vector<size_t> &row_perm, const std::vector<size_t> &column_perm){
	
	using MatT = ColumnMatrix<TC>;
    U_k.permute_rows(bats::util::inv_perm(column_perm));
	auto F = UQL(U_k);                                                
    MatT U = F.U;        
    MatT L = F.L;
    MatT Q = F.E;

    MatT R_k_prime = R_k * l_inv(L) * Q.T();
    R_k_prime.permute_rows(row_perm);
    MatT U_double_prime = MatT::identity(R_k_prime.ncol());
    bats::reduce_matrix(R_k_prime, U_double_prime);

    MatT U_k_prime = U * U_double_prime;

    R_k = R_k_prime;
    U_k = U_k_prime;    
}

// Return an identity permutation, k is its length
std::vector<size_t> identity_perm(const size_t& k){
    std::vector<size_t> l(k);
    std::iota(l.begin(), l.end(), 0);
    return l;
}

/*
Find the index of sorted elements in vector v 
*/
template <typename T>
std::vector<size_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

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
which can be use to determine permutation/addition/deletion.
A: index in the original list of simplices
B: index in the new list of simplices
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
bats::NO_IND is used to denote a simplex not appear in A or B. 
*/
struct Updating_filtration_info_k{
    std::vector<simplex_info> intersection_of_simplices;
    std::vector<size_t> permutation_of_intersection;
    std::vector<simplex_info> deletion_of_simplices;
    std::vector<simplex_info> addition_of_simplices;
    size_t k; // dimension of the simplices

    Updating_filtration_info_k(const std::vector<std::vector<size_t>> &A, 
    const std::vector<std::vector<size_t>> &B, const size_t & dim){
        k = dim;
        // step 1: find the permuation and deletion information
        std::vector<size_t> list_of_index_in_B; // store the indices of itersecting simplcies in B
        for (size_t i = 0; i < A.size(); i++){ // find A \cap B and A - (A \cap B) 
            auto it = std::find(B.begin(), B.end(), A[i]);
            
            if (it != std::end(B)){ // permuation information, since A[i] is finded in B successfully
                auto index_it = std::distance(B.begin(), it); // find the index of finded simplex in B
                // std::cout << "B contains ";
                // print_simplex(*it);
                // std::cout << ", and its index in A is "<< i << ", and in B is " << index_it << std::endl;
                simplex_info inter_splx(*it, i, index_it);
                intersection_of_simplices.emplace_back(inter_splx);
                list_of_index_in_B.emplace_back(index_it);
            }
            else{ // deletion information, since A[i] is not finded in B
                // std::cout << "B does not contain ";
                // print_simplex(A[i]);
                // std::cout << "\n";
                simplex_info delete_splx(A[i], i, bats::NO_IND);
                deletion_of_simplices.emplace_back(delete_splx);
            } 
        }
        permutation_of_intersection = sort_indexes(list_of_index_in_B);

        // step 2: find the addition information, i.e., find the elements in B but not in A \cap B
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
        if (intersection_of_simplices.empty()){
            std::cout << "no intersection" << std::endl;
        }else{
            std::cout << "need intersection" << std::endl;
            for (auto & v: intersection_of_simplices){
                v.print();
            }
            std::cout << "transformed to permuation: ( ";
            for (auto & i: permutation_of_intersection){
                std::cout << i << " ";
            }
            std::cout << ") " << std::endl;
        }
        
        if (deletion_of_simplices.empty()){
            std::cout << "no deletion" << std::endl;
        }else{
            std::cout << "need deletion" << std::endl;
            for (auto & v: deletion_of_simplices){
                v.print();
            }
        }
        
        if (addition_of_simplices.empty()){
            std::cout << "no addition" << std::endl;
        }else{
            std::cout << "need addition" << std::endl;
            for (auto & v: addition_of_simplices){
                v.print();
            }
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
    std::vector<Updating_filtration_info_k> information; 

    Updating_filtration_info(const simplicial_complex& splx1, 
    const simplicial_complex& splx2 ){
        for (size_t i=0; i<splx1.size(); i++){
            information.emplace_back(Updating_filtration_info_k(splx1[i], splx2[i], i));
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
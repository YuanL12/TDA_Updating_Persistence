/*
header file for PLU factorization
*/
#pragma once
#include <bats.hpp>
#include <iostream>
#include <vector>

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
/*
header file for PLU factorization
*/
#pragma once
#include <bats.hpp>
#include <iostream>


void hello_PLU(){
	std::cout << "hello from PLU.hpp" <<std::endl;
}

template <class TC>
void update_reduction(ColumnMatrix<TC> &R_k, ColumnMatrix<TC> &U_k, 
const ColumnMatrix<TC> & row_perm, 
const ColumnMatrix<TC> &column_perm){
	
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

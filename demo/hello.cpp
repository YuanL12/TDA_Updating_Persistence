#include <bats.hpp>
#include <plu.hpp>
#include <vector>

#define F2 ModP<int, 2>
#define F3 ModP<int, 3>

int main() {

    // you can use BATS
    //bats::SimplicialComplex X;
    //X.add_recursive({0,1});
    //X.print_summary();
    
	// How to initialzie a sparse matrix
    //std::vector<size_t> colptr{0, 2, 3, 6};
    //std::vector<size_t> rowind{0, 2, 2, 0, 1, 2};
    //std::vector<int> val{1, 2, 3, 4, 5, 6};
    //CSCMatrix<int> M(3,3,colptr, rowind, val);
    //M.print(0,3,0,3);
    
	
	using VT = SparseVector<F2, size_t>;
	using MatT = ColumnMatrix<VT>;
	//sparse vectos are specified by index and its value(nonzero) 
	std::vector<VT> cols{VT({0},{1}), VT({1,2},{1,1}), VT({0,2},{1,1})};
	
	//intialize a matrix
	MatT B(4, 3, cols);
	//B.swap_cols(0,1);//swap columns
	//B.print();//print columns
	MatT R(B);
	MatT U = MatT::identity(B.ncol());
	bats::reduce_matrix(R, U);

	//check all results
	if (U.is_upper())
		if (R.is_reduced())
			if (R * u_inv(U) == B)
				std::cout << "before permuation, success"<< '\n';

	//swap columns of row permutation matrix will result row change of 
	//the original one
	MatT row_perm(MatT::identity(R.nrow()));
	row_perm.swap_cols(1,2);
	
	//swap columns of column permutation matrix will result same 
	//effect on the column of original matrix
	MatT col_perm(MatT::identity(R.ncol()));
    col_perm.swap_cols(0,1);
	
	update_reduction<VT>(R, U, row_perm, col_perm);

	if (R * u_inv(U) == row_perm * B * col_perm){
		std::cout << "after permuation, success"<< '\n';
	}
	return 0;
}

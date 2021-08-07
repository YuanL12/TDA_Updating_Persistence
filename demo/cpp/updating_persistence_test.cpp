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
    
	using VT = SparseVector<F2, size_t>;
	using MatT = ColumnMatrix<VT>;
	//sparse vectos are specified by index and its value(nonzero) 
	std::vector<VT> cols{VT({0},{1}), VT({1,2},{1,1}), VT({0,2},{1,1})};

	{
		std::cout << "\nstart matrix permuation" << std::endl;
		MatT B(4, 3, cols); //initialize a boundary matrix
		MatT R(B);
		MatT U = MatT::identity(B.ncol());
		bats::reduce_matrix(R, U);

		//check all results
		if (U.is_upper())
			if (R.is_reduced())
				if (R * u_inv(U) == B)
					std::cout << "before permutation, success"<< '\n';

		//swap columns of row permutation matrix will result row change of 
		//the original one
		MatT row_perm_mat(MatT::identity(R.nrow()));
		row_perm_mat.swap_cols(1,2);
		std::vector<size_t> row_perm{0,2,1};
		//swap columns of column permutation matrix will result same 
		//effect on the column of original matrix
		MatT col_perm_mat(MatT::identity(R.ncol()));
		col_perm_mat.swap_cols(0,1);
		std::vector<size_t> col_perm{1,0,2};
		std::cout << "before matrix permutation R is" << std::endl;
		R.print();
	
		update_reduction<VT>(R, U, row_perm_mat, col_perm_mat);
		if (R * u_inv(U) == row_perm_mat * B * col_perm_mat){
			std::cout << "using matrix permutation, success"<< '\n';
		}
		std::cout << "after matrix permutation R is" << std::endl;
		R.print();
	}

	{
		std::cout << "\nstart list permuation" << std::endl;
		MatT B(4, 3, cols); //initialize a boundary matrix
		MatT R(B);
		MatT U = MatT::identity(B.ncol());
		bats::reduce_matrix(R, U);

		//check all results
		if (U.is_upper())
			if (R.is_reduced())
				if (R * u_inv(U) == B)
					std::cout << "before permutation, success"<< '\n';

		//swap columns of row permutation matrix will result row change of 
		//the original one
		std::vector<size_t> row_perm{0,2,1};
		//swap columns of column permutation matrix will result same 
		//effect on the column of original matrix
		std::vector<size_t> col_perm{1,0,2};

		std::cout << "before list permutation R is" << std::endl;
		R.print();
		update_reduction<VT>(R, U, row_perm, col_perm);
		B.permute_cols(col_perm);
		B.permute_rows(row_perm);
		if (R * u_inv(U) == B ){
			std::cout << "using list permutation, success"<< '\n';
		}
		std::cout << "after list permutation R is" << std::endl;
		R.print();
	}

	return 0;
}

#include <bats.hpp>
#include <plu.hpp>
#include <vector>
#include <iostream>
#define F2 ModP<int, 2>

int main()
{
	// How to initialize a sparse matrix
    std::vector<size_t> colptr{0, 2, 3, 6};
    std::vector<size_t> rowind{0, 2, 2, 0, 1, 2};
    std::vector<int> val{1, 2, 3, 4, 5, 6};
    CSCMatrix<int> M(3,3,colptr, rowind, val);
    M.print(0,3,0,3); // print selected parts

    //initialize column matrix by sparse vector in mod 2
    using VT = SparseVector<F2, size_t>;
	using MatT = ColumnMatrix<VT>;
	//sparse vectos are specified by index and its value(nonzero) 
	std::vector<VT> cols{VT({0},{1}), VT({1,2},{1,1}), VT({0,2},{1,1})};
	// creat column
	MatT B(4, 3, cols);
	//B.swap_cols(0,1);//swap columns
	//B.print();//print columns

    MatT C(M); // CSC matrix convert to Column matrix with mod 2
    C.print();
    //std::cout << 1 << std::endl;
    return 0;
}
#include <bats.hpp>
#include <util/print.hpp>
#include <vector>
#include <iostream>
#define F2 ModP<int, 2>

int main()
{
	// Initialize a CSC sparse matrix
    std::vector<size_t> colptr{0, 2, 3, 6};
    std::vector<size_t> rowind{0, 2, 2, 0, 1, 2};
    std::vector<int> val{1, 2, 3, 4, 5, 6};
    CSCMatrix<int> M(3,3,colptr, rowind, val);
    std::cout << "\nMatrix M is" << std::endl;
    M.print(0,3,0,3); // print selected parts

    //Initialize column matrix by sparse vector in mod 2
    using VT = SparseVector<F2, size_t>;
	using MatT = ColumnMatrix<VT>;
	// sparse vectos are specified by index and its value(nonzero) 
    // creat column
	std::vector<VT> cols{VT({0,1,2},{1,1,1}), VT({1,2},{1,1}), VT({0,2},{1,1})};

	MatT B(4, 3, cols);
	//B.swap_cols(0,1);//swap columns
	//B.print();//print matrix

    std::cout << "\nCSC matrix M convert to Column matrix(mod 2) C" << std::endl;
    MatT C(M); // CSC matrix convert to Column matrix with mod 2
    C.print();
    
    // BATs also provide inverse permutation
    std::cout << "\npermutation and its inverse" << std::endl;
    std::vector<size_t> perm = {2,1,0,3};
    auto inv_perm = bats::util::inv_perm(perm);
    std::cout << "original" << std::endl;
    bats::print_1D_vectors(perm);
    std::cout << "its inverse is" << std::endl;
    bats::print_1D_vectors(inv_perm);

    // if you want to permute rows 
    perm = {2,0,1}; // two-line notation, not cycle notation 
    C.permute_rows(perm);
    std::cout << "\npermute the rows by (2 0 1), now C is" << std::endl;
    C.print();

    std::cout << "\npermute the 1st and 3rd columns, now C is" << std::endl;
    perm = {2,1,0};
    C.permute_cols(perm);
    C.print();



    // append a column to the end of a matrix
    C.append_column(VT({0,1,2},{1,1,1}));
    C.print();

    // insert a row a the index 1
    C.insert_column(1, VT({0,1,2},{1,1,1}));
    C.print();

    // append a row to the end of a matrix
    C.append_row({1,3,1,2});
    C.print();

    // insert a row a the position 0
    C.insert_row(0, {0,0,1,1});
    C.print();

    // delete a column with a given index
    std::cout << "\ndelete column" << std::endl;
    C.erase_column(1);
    C.print();

    // delete a row with a given index
    std::cout << "\ndelete row" << std::endl;
    C.erase_row(2);
    C.print();

    std::cout << "\ndelete the first row" << std::endl;
    C.erase_row(0);
    C.print();   

    return 0;
}
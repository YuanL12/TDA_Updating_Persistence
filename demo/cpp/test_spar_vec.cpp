#include <bats.hpp>
#include <plu.hpp>
#include <vector>
#include <iostream>
#define F2 ModP<int, 2>
int main()
{
    using VT = SparseVector<int, size_t>;
    std::vector<size_t>  perm{2,0,1};

	//sparse vectos are specified by index and its value(nonzero) 
    auto v_2 = VT({0,2},{5,6});
    // for (auto& v: v_1.indval){

    // }
    v_2.print();
    v_2.permute(perm);

    std::cout << "\nafter permutation" << std::endl;
    v_2.print();
    
    return 0;
}
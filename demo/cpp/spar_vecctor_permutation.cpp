#include <bats.hpp>
#include <vector>
#include <iostream>
#define F2 ModP<int, 2>
int main()
{
    using VT = SparseVector<int, size_t>;
    std::vector<size_t>  perm{2,0,1};
    std::cout << "\nthe last element of perm is" << std::endl;
    std::cout << *(--perm.end())<< std::endl;
    std::cout << *perm.end() << std::endl;
	//sparse vectos are specified by index and its value(nonzero) 
    std::cout << "\nfor v_1" << std::endl;
    auto v_1 = VT({0,2},{5,6}); // vector [5 0 6]
    v_1.print();
    v_1.permute(perm); // should be [6 5 0]

    std::cout << "after permutation" << std::endl;
    v_1.print();
    if (v_1 == VT({0,1},{6,5})){std::cout << "permutation success" << std::endl;}

    std::cout << "\nfor v_2" << std::endl;
    auto v_2 = VT({0,2,4},{5,6,1}); // vector [5 0 6 0 1]
    v_2.print();
    perm = {2, 0, 4, 3 ,1};
    v_2.permute(perm); // should be [6 5 1 0 0]
    std::cout << "after permutation" << std::endl;
    v_2.print();
    if (v_2 == VT({0,1,2},{6,5,1})){std::cout << "permutation success" << std::endl;}

    {
        perm = {2, 0, 4, 3 ,1};
        auto v_3 = VT({0,2,4},{5,6,1}); 
        std::cout << "\nv_3 is" << std::endl;
        v_3.print();

        for (size_t i = 0; i<perm.size(); i++){
            // find the indval with indval.ind == perm[i] 
            auto it = std::lower_bound(v_3.nzbegin(), v_3.nzend(), perm[i]);
            if (perm[i] == (*it).ind){
                (*it).ind = i; // set new indval.ind = i
            }   
        }
        std::cout << "\nbefore sort" << std::endl;
        v_3.print();
        v_3.sort();
        std::cout << "\nafter sort" << std::endl;
        v_3.print();
    }

    {
        perm = {2, 0, 4, 3 ,1};
        auto v_3 = VT({0,2,4},{5,6,1}); 
        std::cout << "\nv_3 is" << std::endl;
        v_3.print();
        v_3.permute(perm);
        std::cout << "\nafter permutation" << std::endl;
        v_3.print();

    }
    
    // for (auto i = 0; i < 6; i++ ){
    //     auto it = std::lower_bound(v_3.nzbegin(), v_3.nzend(), size_t(i));
    //     if (size_t(i) == (*it).ind)
    //     std::cout << "for i = "<< i<<", *it = "<< *it << ", and (*it).ind = "<< (*it).ind<< std::endl;
    // }    

    return 0;
}
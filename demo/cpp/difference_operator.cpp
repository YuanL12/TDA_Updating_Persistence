#include <iostream>
#include <vector>
int main()
{
    
    std::vector<size_t>  perm{2,0,1};

    for (auto iter = perm.begin(); iter < perm.end(); iter++) std::cout << *iter<<" " << std::endl;
    std::cout << "\n" << std::endl;
    for (auto iter = perm.begin(); iter < perm.end(); ++iter) std::cout << *iter<<" " << std::endl;
    std::cout << "\n" << std::endl;
    for (size_t i =0; i < perm.size(); i++) std::cout << perm[i]<<" " << std::endl;

    
    std::vector<int> v{ 10, 20, 30, 30, 30, 40, 50 };
 
    // Print vector
    std::cout << "Vector contains :";
    for (unsigned int i = 0; i < v.size(); i++)
        std::cout << " " << v[i];
    std::cout << "\n";
 
    std::vector<int>::iterator low1;
     
    // std :: lower_bound
    low1 = std::lower_bound(v.begin(), v.end(), 70);
    std::cout << *low1 << std::endl;
    return 0;
}
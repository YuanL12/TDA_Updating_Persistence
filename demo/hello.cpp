#include <bats.hpp>
#include <plu.hpp>

int main() {

    // you can use BATS
    bats::SimplicialComplex X;
    X.add_recursive({0,1});
    X.print_summary();

    // you can use functions from your
    // repository include headers
    hello_PLU();

    return 0;
}

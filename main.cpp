#include <chrono>
#include <iostream>
#include "graph/tree/isomorphism.h"
#include "nt/dirichelet.h"
#include <cstdint>


int main()
{
    using cp::integer;
    integer n;
    std::cin >> n;
    std::vector<integer> H(5);
    for (integer i = 0; i < 5; i++)
        H[i] = cp::pow(3LL, i*n);
    std::vector<integer> µ = {1,-1};
    auto J = cp::convolve(H,µ);
    for (auto j:J) std::cout << j << ' ';
}

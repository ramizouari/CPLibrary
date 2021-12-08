//
// Created by ramizouari on 01/12/2021.
//
#include "linear_algebra.h"
#include <iostream>
#include "abstract_algebra.h"
#include "fft.h"
#include "modular_arithmetic.h"
#include "optimisation.h"

int main()
{
    d_cyclic::m=1e9+9;
    std::vector<d_cyclic> u={1,1,1,1,1,1,1,1,1,1,1};
    factoriser F(2e5);
    for(auto s: fast_multiply(u,u,F))
        std::cout << (integer)s << ' ';
}
//
// Created by ramizouari on 01/12/2021.
//
#include "linear_algebra.h"
#include "modular_arithmetic.h"
#include <iostream>
#include "fft.h"
#include "ring_extension.h"


int main()
{
    s_matrix<real,3,3> M({{3,0,2},{0,0,0},{5,0,1}});
    s_vector<real,3> v({21,0,7});
    for(auto s:M.solve(v))
        std::cout << s << ' ';
}
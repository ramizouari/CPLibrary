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
    ring_extension<real>::q=polynomial<real>({1,1,1});
    ring_extension<real> p({0,1});
    for(auto s:pow(p,1001))
        std::cout << s << ' ';
}
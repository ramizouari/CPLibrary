//
// Created by ramizouari on 01/12/2021.
//
#include "linear_algebra.h"
#include <iostream>
#include "abstract_algebra.h"
#include "fft.h"
#include "ring_extension.h"
#include <chrono>
constexpr integer m=1e9+7;
using IK=cyclic<m>;
using IQ=rational_extension<integer>;
int main()
{
    polynomial<IQ> p=newton_interpolation<IQ>({1,2,3,4,5,6},{3,9,30,101,358,1443});
    integer n;
    std::cin >> n;
    rational_t<integer> A=p(IQ(n));
    auto [u,v]=A;
    std::cout << u << ' ' << v;
}
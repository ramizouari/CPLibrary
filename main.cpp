//
// Created by ramizouari on 01/12/2021.
//
#include "linear_algebra.h"
#include <iostream>
#include "abstract_algebra.h"
#include "fft.h"
#include "ring_extension.h"
#include <chrono>
#include "optimisation.h"
using E=d_vector<real>;

namespace global
{
    d_matrix<real> M({{5,1},{1,2}});
    d_vector<real> b({3,5});
}

real f(E a)
{
    using namespace global;
    static L2_inner_product<E> B;
    return B.inner_product(a,M*a+b);
}

int main()
{
    derivator<E,real,E> D;
    barzilai_borwein_gradient_descent<E> GD(E({0,0}),D,.1);
    auto u=GD.argmin(f);
    for(auto s:u)
        std::cout << s << ' ';
    using namespace global;
    std::cout << '\n';
    for(auto s:(-.5L*M.inv())*b)
        std::cout << s << ' ';
}
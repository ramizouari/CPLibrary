//
// Created by ramizouari on 15/04/22.
//
#include <iostream>
#include "nt/modular_arithmetic.h"
#include "polynomial/polynomial.h"
constexpr integer M=998'244'353;
using IF=cyclic<M>;

int main()
{
    int n,&N=n;
    std::cin >> n;
    std::vector<IF> A(n),B(n);
    for(auto &a:A)
        std::cin >> (integer&)a;
    for(auto &b:B)
        std::cin >> (integer&)b;
    polynomial<IF> P(A),Q(B);
    std::vector<IF> H(N),I(N);
    for(int i=0;i<N;i++)
    {
        I[i]=i;
        H[i]=P(Q(IF{i}));
    }
    auto C=newton_interpolation(I,H);
    for(auto c:C)
        std::cout << (integer)c << ' ';
}
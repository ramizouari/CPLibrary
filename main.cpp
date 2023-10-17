#include <iostream>
#include <vector>
#include "nt/number_theory.h"
#include "nt/modular_arithmetic.h"
#include "nt/primality.h"
#include "boost/multiprecision/cpp_int.hpp"

using big_integer = boost::multiprecision::cpp_int;

constexpr integer L=2e6, M=998244353;
using IK=cyclic<M>;


std::vector<std::pair<integer,integer>> prime_decomposition(integer x)
{
    auto L=std::sqrt(x);
    std::vector<std::pair<integer,integer>> D;
    for(int i=2;i<=L;i++)if (x%i==0)
    {
        D.emplace_back(i,0);
        while(x%i==0)
        {
            D.back().second++;
            x/=i;
        }
    }
    if(x>1)
        D.emplace_back(x,1);
    return D;
}

int main()
{
    factoriser F(L);
    fast_factoriser FF(20,polynomial<integer>({1,0,1}),25);
    randomized_fast_factoriser RFF(20,5);
    integer n,m;
    std::cin >> n >> m;
    auto D=RFF.prime_decomposition(n);
    m%=M;
    IK d=m;
    bool square_free=false;
    for(auto [_,r]:D)
    {
        d *= m * r + 1;
        if(r%2 == 1)
            square_free=true;
    }
    if(m%2 == 1 && !square_free)
        d--;
    d/=2;
    big_integer x=85;
    std::cout << static_cast<integer>(d)+x << '\n';
}
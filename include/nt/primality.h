//
// Created by ramizouari on 16/10/22.
//

#ifndef CPLIBRARY_PRIMALITY_H
#define CPLIBRARY_PRIMALITY_H

#include "algebra/abstract_algebra.h"
#include "modular_arithmetic.h"

bool rabin_miller_primality_test(integer n, integer _a)
{
    if (n <= 2)
        return n == 2;
    else if (n % 2 == 0)
        return false;
    integer r = n - 1, h = 0;
    while (r % 2 == 0)
    {
        r /= 2;
        h++;
    }
    integer d = 1;
    d_cyclic a(_a);
    auto u = pow(a, r);
    if (u - 1 == 0)
        return true;
    for (int i = 0; i <= h; i++, u *= u) if (u + 1 == 0)
            return true;
    return false;
}

bool rabin_miller(integer n, integer iter = 7)
{
    static std::random_device dev;
    static std::mt19937_64 g(dev());
    if (n == 1)
        return false;
    std::uniform_int_distribution<integer> d(2, n - 1);
    d_cyclic::m = n;
    for (int i = 0; i < iter; i++) if (!rabin_miller_primality_test(n, d(g)))
            return false;
    return true;
}

bool rabin_miller(integer n, const std::vector<integer>& provers)
{
    if (n == 1)
        return false;
    d_cyclic::m = n;
    for (const auto& d : provers) if (!rabin_miller_primality_test(n, d))
            return false;
    return true;
}

bool fermat_test(integer n, const std::vector<integer>& provers)
{
    if (n == 1)
        return false;
    d_cyclic::m = n;
    for (const auto &p : provers) if (pow<d_cyclic>(p, n - 1) != 1)
            return false;
    return true;
}

integer rho_divisor_method(integer n,const polynomial<integer> &_P,integer x0)
{
    d_cyclic::m=n;
    polynomial<d_cyclic> P;
    for(int i=0;i<=_P.degree();i++)
        P.p.emplace_back(_P[i]);
    P.reduce();
    d_cyclic x=x0,y=x;
    integer d=1;
    do
    {
        x=P(x);
        y=P(P(y));
        d=gcd(static_cast<integer>(y-x),n);
    }while(d==1);
    return d;
};


class fast_factoriser
{
    int iters;
    polynomial<integer> P;
    integer x0;
public:
    fast_factoriser(int iters,const polynomial<integer> &P,integer x0):iters(iters),P(P),x0(x0){}
    [[nodiscard]] integer smallest_divisor(integer n) const
    {
        if(rabin_miller(n,iters))
            return n;
        return smallest_divisor(rho_divisor_method(n,P,x0));
    }
};

class randomized_fast_factoriser
{
    int iters;
    std::random_device dev;
    std::mt19937_64 g{dev()};
    int max_degree;
public:
    randomized_fast_factoriser(int iters,int max_degree):iters(iters),max_degree(max_degree){}
    [[nodiscard]] integer smallest_divisor(integer n)
    {
        std::uniform_int_distribution<integer> d(0,n-1);
        std::vector<integer> P(max_degree+1);
        std::generate(P.begin(),P.end(),std::bind(d,g));
        integer x0=d(g);
        if(rabin_miller(n,iters))
            return n;
        return smallest_divisor(rho_divisor_method(n,P,x0));
    }
};


#endif //CPLIBRARY_PRIMALITY_H

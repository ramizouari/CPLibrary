//
// Created by ramizouari on 16/10/22.
//

#ifndef CPLIBRARY_PRIMALITY_H
#define CPLIBRARY_PRIMALITY_H

#include "algebra/abstract_algebra.h"
#include "modular_arithmetic.h"

inline bool rabin_miller_primality_test(integer n, integer _a)
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


#endif //CPLIBRARY_PRIMALITY_H

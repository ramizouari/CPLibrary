//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_FUNCTIONS_H
#define CPLIBRARY_FUNCTIONS_H
#include "number_theory.h"
namespace cp
{
    inline integer totient(integer n,abstract_factoriser &F)
    {
        integer r=1;
        for(auto [p,k]:F.prime_decomposition(n))
            r*=pow(p,k-1)*(p-1);
        return r;
    }

    inline integer carmichael_totient(integer n,abstract_factoriser &F)
    {
        integer r=1;
        for(auto [p,k]:F.prime_decomposition(n))
        {
            if(p==2&&k>=3)
                r=std::lcm(r,pow(2,k-2));
            else
                r = std::lcm(r, pow(p, k - 1) * (p - 1));
        }
        return r;
    }

    inline integer divisors_count(integer n,abstract_factoriser &F)
    {
        integer r=1;
        for(auto [p,k]:F.prime_decomposition(n))
            r*=k+1;
        return r;
    }

    inline integer divisors_sum(integer n,abstract_factoriser &F)
    {
        integer r=1;
        for(auto [p,k]:F.prime_decomposition(n))
            r*=(pow(p,k+1)-1)/(p-1);
        return r;
    }

    inline integer prime_multiplicity(integer n,integer p,abstract_factoriser &F)
    {
        for(auto [q,k]:F.prime_decomposition(n))
            if(q==p)
                return k;
        return 0;
    }

    inline integer prime_multiplicity(integer n,integer p)
    {
        integer r=0;
        while(n%p==0)
        {
            n/=p;
            r++;
        }
        return r;
    }


}

#endif //CPLIBRARY_FUNCTIONS_H

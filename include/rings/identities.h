//
// Created by ramizouari on 15/12/24.
//

#ifndef CPLIBRARY_RING_IDENTITIES_H
#define CPLIBRARY_RING_IDENTITIES_H

#include "algebra/abstract_algebra.h"
#include "algebra/structures.h"
#include <concepts>

namespace cp::rings
{
    /*
     * Compute sum x^k, a <= k < b
     */
    template<ring R>
    R geometric_power_sum(const R& x, integer a,integer b)
    {
        if (a) return pow(x,a) * geometric_power_sum(x,0,b-a);
        if (b==0) return {};
        auto c=(a+b)/2;
        auto w = pow(x,c);
        auto y1 = geometric_power_sum(x,0,c);
        auto y2 = y1;
        if (b % 2) y2 += w;
        return y1 + w * y2;
    }

    template<ring R, ring  integer_like= R> requires (std::convertible_to<integer_like,R>)
    std::vector<R> geometric_power_sum_rec(const R& x, integer n, std::size_t r,const std::vector<std::vector<integer_like>> & nCr)
    {
        if (n==0) return std::vector<R>(r+1);
        auto m =n/2;
        auto w = pow(x,m);
        integer_like k=1,s=1;
        auto psi = geometric_power_sum_rec(x,m,r,nCr);
        std::vector<R> phi(r+1);
        for (int p=0;p<=r;p++) {
            phi[p] += psi[p];
            if (n%2) psi[p]+= k*w;
            for (int i=0;i<=p;i++) {
                phi[p] +=  s*nCr[p][p-i] * w * psi[p-i];
                s*=m;
            }
            k*=m;
            s=1;
        }
        return phi;
    }


    /*
     * Compute sum k^r*x^k, a <= k < b,
     */
    template<ring R, ring  integer_like= R> requires (std::convertible_to<integer_like,R>)
    R geometric_power_sum(R x, integer a,integer b, std::size_t r)
    {
        std::vector<std::vector<integer_like>> nCr(r+1,std::vector<integer_like>(r+1));
        nCr[0][0]=1;
        for (int i=0;i<=r;i++) for (int j=0;j<=r;j++) {
            if (i) nCr[i][j]+=nCr[i-1][j];
            if (i && j) nCr[i][j]+=nCr[i-1][j-1];
        }
        return geometric_power_sum_rec(x,b,r,nCr).back() - geometric_power_sum_rec(x,a,r,nCr).back();
    }
};

#endif //CPLIBRARY_RING_IDENTITIES_H

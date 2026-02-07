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

    template<ring R>
    void build_powers(const R &x, long long n, std::vector<R> &W)
    {
        if (n==0) {
            W.emplace_back(1);
            return;
        }
        build_powers(x,n/2,W);
        auto w = W.back();
        w*=w;
        if (n%2) w*=x;
        W.push_back(w);
    }

    /*
     * Compute sum x^k, a <= k < b
     */
    template<ring R>
    R geometric_power_sum(const R& x, integer n, std::vector<R> &W, unsigned depth=0)
    {
        if (depth==W.size()) return {};
        auto m=n/2;
        auto w = W[W.size()-1-depth];
        auto y1 = geometric_power_sum(x,m,W,depth+1);
        auto y2 = y1;
        if (n % 2) y2 += w;
        return y1 + w * y2;
    }

    template<ring R>
    R geometric_power_sum(const R& x, integer a,integer b)
    {
        if (a>=b) return {};
        std::vector<R> W;
        build_powers(x,(b-a)/2,W);
        if (a) return pow(x,a) * geometric_power_sum(x,b-a,W);
        return geometric_power_sum(x,b-a,W);
    }

    template<ring R, ring  integer_like= R> requires (std::convertible_to<integer_like,R>)
    std::vector<R> geometric_power_sum_rec(const R& x, integer n, std::size_t r,const std::vector<std::vector<integer_like>> & nCr,
        const std::vector<R> & W, unsigned depth=0)
    {
        if (depth==W.size()) return std::vector<R>(r+1);
        auto m =n/2;
        auto w = W[W.size()-depth-1];
        integer_like k=1,s=1;
        auto psi = geometric_power_sum_rec(x,m,r,nCr,W,depth+1);
        std::vector<R> phi(r+1);
        for (int p=0;p<=r;p++) {
            phi[p] += psi[p];
            if (n%2) psi[p]+= k*w;
            for (int i=0;i<=p;i++) {
                phi[p] +=  s*nCr[p][p-i] * w * psi[p-i];
                s*=integer_like(m);
            }
            k*=integer_like(m);
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
        std::vector<R> Wa,Wb;
        build_powers(x,b/2,Wb);
        build_powers(x,a/2,Wa);
        return geometric_power_sum_rec(x,b,r,nCr,Wb).back() - geometric_power_sum_rec(x,a,r,nCr,Wa).back();
    }
};

#endif //CPLIBRARY_RING_IDENTITIES_H

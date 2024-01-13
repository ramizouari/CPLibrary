//
// Created by ramizouari on 12/01/24.
//

#ifndef CPLIBRARY_FORMAL_SERIES_H
#define CPLIBRARY_FORMAL_SERIES_H
#include "fast_polynomial.h"

namespace cp
{
    // Formal Newton method. Requires z0 to be a solution of P(z)=0 in R[x]/x
    template<typename R,typename X,typename Functional,typename Derivative>
    polynomial<R> formal_newton_method_2(polynomial<R> P, const X &z0, int n,const Functional& F, const Derivative & dF)
    {
        if(n==1)
            return z0;
        else
        {
            auto z=formal_newton_method_2(P,z0,n/2,F,dF);
            auto w=F(z,n),dw_inv=formal_inv(dF(z,n),n);
            auto dz=fast_multiplication(w, dw_inv);
            dz.data().resize(n);
            return z - dz;
        }
    }

    template<typename R,typename X,typename Functional,typename Derivative>
    polynomial<R> formal_newton_method(polynomial<R> P, const X &z0, int n,const Functional& F, const Derivative & dF)
    {
        auto Q= formal_newton_method_2(P,z0,std::bit_ceil<unsigned>(n),F,dF);
        Q.data().resize(n);
        return Q;
    }

    // Formal logarithm
    template<typename R>
    polynomial<R> formal_logarithm(const polynomial<R> &P,int n)
    {
        auto Q= fast_multiplication(P.derivative(),formal_inv(P,n));
        for(int i=Q.data().size();i>0;i--)
        {
            Q[i] = Q[i - 1];
            Q[i] /= i;
        }
        if(!Q.data().empty())
            Q[0]=0;
        Q.data().resize(n);
        return Q;
    }

    template<typename R>
    std::vector<R> formal_logarithm(const std::vector<R> &P,int n)
    {
        auto dP=P;
        dP.erase(dP.begin());
        for(int i=0;i<dP.size();i++)
            dP[i]*=i+1;
        auto Q= fast_multiplication(dP,formal_inv(P,n));
        Q.resize(n);
        for(int i=Q.size();i>0;i--)
        {
            Q[i] = Q[i - 1];
            Q[i] /= i;
        }
        if(!Q.empty())
            Q[0]=0;
        return Q;
    }

    template<typename R>
    polynomial<R> formal_exp_2(const polynomial<R> &P,int n)
    {
        if(n==1)
            return 1;
        else
        {
            auto Q=formal_exp_2(P,n/2);
            Q=fast_multiplication(Q,P + R(1) - formal_logarithm(Q,n));
            Q.data().resize(n);
            return Q;
        }
    }

    template<typename R>
    polynomial<R> formal_exp(const polynomial<R> &P,int n)
    {
        auto Q= formal_exp_2(P,std::bit_ceil<unsigned>(n));
        Q.data().resize(n);
        return Q;
    }

}

#endif //CPLIBRARY_FORMAL_SERIES_H

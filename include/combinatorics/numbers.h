//
// Created by ramizouari on 11/01/24.
//

#ifndef CPLIBRARY_COMBINATORICS_FUNCTIONS_H
#define CPLIBRARY_COMBINATORICS_FUNCTIONS_H
#include "factorial.h"
#include "binomial.h"
namespace cp
{


    template<typename T>
    T nCr(integer n,integer r,const abstract_nCr<T> &C)
    {
        return C(n,r);
    }

    template<typename T>
    T nCr(integer n,integer m)
    {
        return nCr(n,m,naive_nCr<T>());
    }

    template<typename T>
    T nPr(integer n,integer r,const abstract_nCr<T> &C,const abstract_factorial<T> &F)
    {
        return C(n,r)*F(r);
    }

    template<typename T>
    T nPr(integer n,integer r,const abstract_factorial<T> &F,const abstract_inverse_factorial<T> &F_inv)
    {
        return F(n)*F_inv(n-r);
    }

    template<typename T>
    T nPr(integer n,integer r,const abstract_factorial<T> &F)
    {
        return F(n)/F(n-r);
    }

    template<typename T>
    T nPr(integer n,integer r)
    {
        return nPr(n,r,naive_factorial<T>());
    }

    template<typename T>
    T catalan_number(integer n,const abstract_nCr<T> &C)
    {
        return nCr(2*n,n,C)/(n+1);
    }

    template<typename T>
    T catalan_number(integer n, const abstract_nCr<T> &C, const std::vector<T> &I)
    {
        return nCr(2*n,n,C)*I[n];
    }

    template<typename T>
    T catalan_number(integer n)
    {
        return catalan_number(n,naive_nCr<T>());
    }

    template<typename T>
    T multinomial(const std::vector<integer> &A,const abstract_factorial<T> &F)
    {
        integer n=std::accumulate(A.begin(),A.end(),0);
        T res=F(n);
        for(auto a:A)
            res/=F(a);
        return res;
    }

    template<typename T>
    T multinomial(const std::vector<integer> &A,const abstract_factorial<T> &F,const abstract_inverse_factorial<T> &F_inv)
    {
        integer n=std::accumulate(A.begin(),A.end(),0);
        T res=F(n);
        for(auto a:A)
            res*=F_inv(a);
        return res;
    }

    template<typename T>
    T multiset_coefficient(integer n,integer r,const abstract_nCr<T> &C)
    {
        return C(n+r-1,r);
    }

    template<typename T>
    T multiset_coefficient(integer n,integer r)
    {
        return multiset_coefficient(n,r,naive_nCr<T>());
    }
}
#endif //CPLIBRARY_FUNCTIONS_H

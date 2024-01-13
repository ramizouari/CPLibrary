//
// Created by ramizouari on 12/01/24.
//

#ifndef CPLIBRARY_STIRLING_H
#define CPLIBRARY_STIRLING_H

#include "algebra/abstract_algebra.h"

namespace cp
{
    template<typename T>
    struct abstract_stirling_1
    {
        virtual ~abstract_stirling_1()= default;
        virtual T stirling_1(integer n,integer k) const=0;
        T operator()(integer n,integer k) const
        {
            return stirling_1(n,k);
        }
    };

    template<typename T>
    struct abstract_stirling_2
    {
        virtual ~abstract_stirling_2()= default;
        virtual T stirling_2(integer n,integer k) const=0;
        T operator()(integer n,integer k) const
        {
            return stirling_2(n,k);
        }
    };

    template<typename T>
    struct precomputed_stirling_1 : public abstract_stirling_1<T>
    {
        std::vector<std::vector<T>> S;
        explicit precomputed_stirling_1(integer n)
        {
            S.resize(n+1);
            for(int i=0;i<=n;i++)
            {
                S[i].resize(i+1);
                S[i][0]=0;
                S[i][i]=1;
            }
            for(int i=1;i<=n;i++) for(int j=1;j<i;j++)
                S[i][j]=S[i-1][j-1]+(i-1)*S[i-1][j];
        }
        T stirling_1(integer n,integer k) const override
        {
            return S[n][k];
        }
    };

    template<typename T>
    struct precomputed_stirling_2 : public abstract_stirling_2<T>
    {
        std::vector<std::vector<T>> S;
        explicit precomputed_stirling_2(integer n)
        {
            S.resize(n+1);
            for(int i=0;i<=n;i++)
            {
                S[i].resize(i+1);
                S[i][0]=0;
                S[i][i]=1;
            }
            for(int i=1;i<=n;i++) for(int j=1;j<i;j++)
                S[i][j]=S[i-1][j-1]+j*S[i-1][j];
        }
        T stirling_2(integer n,integer k) const override
        {
            return S[n][k];
        }
    };
}
#endif //CPLIBRARY_STIRLING_H

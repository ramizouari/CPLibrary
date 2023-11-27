//
// Created by ASUS on 01/12/2021.
//

#ifndef ACPC_PREPARATION_FAST_POLY_H
#define ACPC_PREPARATION_FAST_POLY_H
#include "fft.h"
#include "ring_extension.h"
#include "algebra/binary_operation.h"
#include "data_structures/range_queries.h"

namespace cp
{
    template<typename R>
    std::vector<R> formal_inv_2(const std::vector<R> &A,int m)
    {
        if(m==1)
            return {R(1)/A.front()};
        auto B=A;
        for(int i=1;i<A.size();i+=2)
            B[i]=-B[i];
        auto C= fast_multiplication(A,B);
        std::vector<R> T;
        T.resize(m/2);
        for(int i=0;i<T.size() && 2*i < C.size();i++)
            T[i]=C[2*i];
        auto S=formal_inv_2(T,m/2);
        std::vector<R> Q;
        Q.resize(m);
        for(int i=0;i<m/2;i++)
            Q[2*i]=S[i];
        return fast_multiplication(B, Q);
    }

    template<typename R>
    std::vector<R> formal_inv(const std::vector<R> &A,int m)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m));
        C.resize(m);
        return C;
    }

    template<typename R>
    polynomial<R> formal_inv(const polynomial<R> &A,int m)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m));
        C.resize(m);
        return C;
    }

    template<typename R,int nilpotence>
    nilpotent_extension<R,nilpotence> fast_inv(const nilpotent_extension<R,nilpotence> &A)
    {
        return nilpotent_extension<R,nilpotence>(formal_inv(A.p,nilpotence));
    }

    template<typename R>
    d_nilpotent_extension<R> fast_inv(const d_nilpotent_extension<R> &A)
    {
        return d_nilpotent_extension<R>(formal_inv(A.p),nilpotence_t{A.nilpotence});
    }


    template<typename R>
    std::vector<R> fast_division(std::vector<R> A,std::vector<R> Q)
    {
        if(A.size()<Q.size())
            return {};
        int m=A.size()-Q.size()+1;
        std::reverse(A.begin(),A.end());
        std::reverse(Q.begin(),Q.end());
        auto P= fast_multiplication(A, formal_inv(Q,m));
        P.resize(m);
        std::reverse(P.begin(),P.end());
        return P;
    }

    template<typename R>
    polynomial<R> fast_division(const polynomial<R> &A,const polynomial<R> &B)
    {
        return fast_division((const std::vector<R>&)A,(const std::vector<R>&)B);
    }

    template<typename R>
    polynomial<R> fast_mod(const polynomial<R>&A,const polynomial<R>& B)
    {
        auto P= fast_division(A,B);
        return A- fast_multiplication(B,P);
    }

    template<typename R>
    std::pair<polynomial<R>,polynomial<R>> fast_euclidean_division(const polynomial<R> &A,const polynomial<R>& B)
    {
        auto P= fast_division(A,B);
        return std::make_pair(P,A- fast_multiplication(B,P));
    }

    template<typename T>
    struct fast_multiplies_t : public multiplies_t<T>
    {
        T reduce(const T&a,const T&b) const override
        {
            return fast_multiplication(a,b);
        }
    };

    template<typename R>
    polynomial<R> fast_polynomial_expansion(const std::vector<R> &X)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        fast_multiplies_t<polynomial<R>>::neutral=1;
        segment_tree<polynomial<R>,fast_multiplies_t<polynomial<R>>> S(P);
        return S.S[0][0];
    }


    template<typename R>
    std::vector<R> fast_multi_evaluation(const polynomial<R> &A,const std::vector<R> &X)
    {
        fast_multiplies_t<polynomial<R>>::neutral=1;
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        segment_tree<polynomial<R>,fast_multiplies_t<polynomial<R>>> S(P);
        std::vector<std::vector<polynomial<R>>> Z(S.h);
        for(int i=0;i<S.h;i++)
            Z[i].resize(1<<i);
        Z[0][0]=fast_mod(A,S.S[0][0]);
        for(int i=1;i<S.h;i++)
            for(int j=0;j<Z[i].size();j++)
                Z[i][j]=fast_mod(Z[i-1][j>>1],S.S[i][j]);
        std::vector<R> Y;
        Y.reserve(n);
        for(int i=0;i<n;i++)
            Y.push_back(Z.back()[i](R{}));
        return Y;
    }
}
#endif //ACPC_PREPARATION_FFT_H

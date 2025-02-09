//
// Created by ramizouari on 15/09/24.
//

#ifndef CPLIBRARY_DIRICHELET_H
#define CPLIBRARY_DIRICHELET_H
#include "algebra/structures.h"
#include "number_theory.h"

namespace cp
{

    template<ring I>
    std::vector<I> conv(std::span<const I> X,std::span<const I> Y)
    {
        auto a=X.size(),b=Y.size(), c =std::max(X.size(),Y.size());
        std::vector<I> Z(c);
        for(int i=0;i<a;i++)for(int j=0;i+j < c && j < b;j++)
            Z[i+j]+=X[i]*Y[j];
        return Z;
    }

    template<ring I>
    std::vector<I> conv_pow(std::span<I> X,unsigned long long n)
    {
        return functional_pow<std::vector<I>>(X,n,[](const auto &x,const auto &y)
        {
            return conv(x,y);
        },std::vector<I>(1,1));
    }


    template<ring I>
    std::vector<I> conv_inverse(std::span<const I> X)
    {
        auto r=X.size();
        std::vector<I> Y(r);
        Y[0]=1;
        for(int i=1;i<r;i++)
        {
            I z{};
            for(int j=0;j < i;j++)
                z+=X[i-j]*Y[j];
            Y[i]=-z;
        }
        return Y;
    }


    template<ring I>
    class convolution_t
    {
        std::span<const I> A,B;
    public:
        convolution_t(std::span<const I> A,std::span<const I> B): A(A),B(B){}
        I operator()(int m,const abstract_factoriser & F)
        {
            auto p = F.smallest_divisor(m);
            I result;
            for(integer q =1; q <= m; q*=p)
                result+= A[m/q] * B[q];
            return result;
        }

        std::vector<I> operator()(const std::vector<int> & J,const abstract_factoriser & F)
        {
            std::vector<I> H(J.size(),1);
            for(int k=1;k<J.size();k++)
                H[k] = (*this)(J[k],F);
            return H;
        }
    };

    template<typename I>
    convolution_t(std::vector<I> & ,std::vector<I> & ) -> convolution_t<I>;
    template<typename I,std::size_t N>
    convolution_t(std::array<I, N> & ,std::array<I,N> &) -> convolution_t<I>;


    template<ring I>
    std::vector<I> convolve(const std::vector<I> & A, const std::vector<I> &B)
    {
        if(A == std::vector<I>{1} )
            return B;
        if(B==  std::vector<I>{1})
            return A;
        auto r= std::max(A.size(),B.size());
        std::vector<I> C(r);
        for(int i=0;i<A.size();i++) for(int j=0;i+j < r && j<B.size() ;j++)
            C[i+j] += A[i]*B[j];
        return C;
    }


    template<ring I>
    class convolution_pow_t
    {
        std::span<const I> X;

        unsigned long long n;

    public:
        convolution_pow_t(std::span<const I> X, unsigned long long n): X(X),n(n)
        {
        }

        std::vector<I> operator()(const std::vector<int> &J, const abstract_factoriser & F)
        {
            std::vector<I> H(J.size());
            for(int k=0;k<J.size();k++) H[k] = X[J[k]];
            return functional_pow(H,n, &convolve<I>,std::vector<I>{1});
        }

    };

    template<ring I>
    convolution_pow_t(std::vector<I> & , unsigned long long) -> convolution_pow_t<I>;
    template<ring I,std::size_t N>
    convolution_pow_t(std::array<I, N> &, unsigned long long) -> convolution_pow_t<I>;


    template<ring I>
    class convolution_inverse_t
    {
        std::span<const I> X;

    public:
        explicit convolution_inverse_t(std::span<const I> X): X(X)
        {
        }

        std::vector<I> operator()(const std::vector<int> &J)
        {
            auto r = J.size();
            std::vector<I> Y(r);
            Y[0]=1;
            for(int i=1;i<r;i++)
            {
                I z{};
                for(int j=0;j < i;j++)
                    z+=X[J[i-j]]*Y[j];
                Y[i]=-z;
            }
            return Y;
        }
    };

    template<ring I>
    convolution_inverse_t(std::vector<I> &) -> convolution_inverse_t<I>;
    template<ring I,std::size_t N>
    convolution_inverse_t(std::array<I, N> &) -> convolution_inverse_t<I>;
}

#endif //CPLIBRARY_DIRICHELET_H

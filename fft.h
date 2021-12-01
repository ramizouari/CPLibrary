//
// Created by ASUS on 01/12/2021.
//

#ifndef ACPC_PREPARATION_FFT_H
#define ACPC_PREPARATION_FFT_H
#include "abstract_algebra.h"
#include <numbers>
#include "polynomial.h"
#include "number_theory.h"
#include <algorithm>

template<typename R>
std::vector<R> apply_powers(const std::vector<R> &A,R w)
{
    R x=1;
    std::vector<R> B(A);
    for(auto &b:B)
    {
        b*=x;
        x*=w;
    }
    return B;
}

template<typename R>
std::vector<R> operator*(R k,const std::vector<R>&T)
{
    std::vector<R> S(T);
    for(auto &s:S)
        s*=k;
    return S;
}

template<typename R>
std::vector<R> operator*(const std::vector<R>&T,R k)
{
    std::vector<R> S(T);
    for(auto &s:S)
        s*=k;
    return S;
}

template<typename R>
std::vector<R> operator/(const std::vector<R>&T,R k)
{
    std::vector<R> S(T);
    for(auto &s:S)
        s/=k;
    return S;
}

template<typename R>
std::vector<R> operator+(const std::vector<R> &A,const std::vector<R> &B)
{
    if(A.size()<B.size())
        return B+A;
    std::vector<R> C(A);
    int m=B.size();
    for(int i=0;i<m;i++)
        C[i]+=B[i];
    return C;
}

template<typename R>
std::vector<R>& operator+=(std::vector<R> &A,const std::vector<R> &B)
{
    if(A.size()<B.size())
        A.resize(B.size());
    int m=B.size();
    for(int i=0;i<m;i++)
        A[i]+=B[i];
    return A;
}


const real pi = acos(-1);

template<bool is_inverse=false>
struct fast_fourier
{
    inline static const factoriser F=factoriser(1e5);
    int n;
    std::complex<real> w;
    using IC=std::complex<real>;
    inline static constexpr int sign=is_inverse?1:-1;
public:
    fast_fourier(int _n):n(_n),w(std::exp(IC(0,2*sign*pi/n)))
    {
    }
    std::vector<IC> unnormalized(const std::vector<IC> &X) const
    {
        if(n==1)
            return X;
        auto p=F.smallest_divisor(n),q=n/p;
        fast_fourier<is_inverse> FFT(q);
        std::vector<std::vector<IC>> U(p,std::vector<IC>(q));
        for(int i=0;i<n;i++)
            U[i%p][i/p]=X[i];
        std::vector<std::vector<IC>> V(p);
        for(int i=0;i<p;i++)
            V[i]=FFT.unnormalized(U[i]);
        std::vector<std::vector<IC>> Q(p,std::vector<IC>(q,0));
        IC z=std::pow(w,q);
        for(int i=0;i<p;i++) for(int j=0;j<p;j++)
                Q[i]+=std::pow(z,i*j)* apply_powers(V[j],std::pow(w,j));
        std::vector<IC> R(n);
        for(int i=0;i<p;i++) for(int j=0;j<q;j++)
                R[i*q+j]=Q[i][j];
        return R;
    }
    std::vector<IC> operator()(const std::vector<IC> &X) const
    {
        return unnormalized(X);
    }
    std::vector<IC> normalized(const std::vector<IC>&X) const
    {
        return  unnormalized(X)/std::sqrt<real>(n);
    }
};

template<typename real>
using inverse_fast_fourier=fast_fourier<true>;


template<int n,typename T>
struct tensor_t
{
    using tensor=std::vector<typename tensor_t<n-1,T>::tensor>;
    std::vector<tensor_t<n-1,T>> U;
    T operator[](const std::array<T,n> &I) const
    {
        std::array<T,n-1> subI;
        for(int i=1;i<n;i++)
            subI[i-1]=I[i];
        return U[subI];
    }
    explicit operator std::vector<tensor_t<n-1,T>>&() const
    {
        return U;
    }
};
template<typename T>
struct tensor_t<0,T>
{
    using tensor=T;
    tensor U;
    T operator[](const std::array<T,0>&)
    {
        return U;
    }
    operator const T&() const
    {
        return U;
    }
};
template<int n,typename T>
using tensor=typename tensor_t<n,T>::tensor;

template<typename T,int n>
T get(const tensor<n,T> &A,std::array<int,n> I)
{
    if constexpr (n==0)
        return A;
    else
    {
        std::array<int, n - 1> subI;
        for(int i=1;i<n;i++)
            subI[i-1]=I[i];
        return get<T,n-1>(A[I[0]],subI);
    }
}
template<typename T,int n>
tensor<n,T> reshape(const std::vector<T> &A,std::array<int,n> shape)
{
    if constexpr (n==0)
        return A[0];
    else
    {
        int m=A.size()/shape[0];
        std::vector<std::vector<T>> B(shape[0],std::vector<T>(m));
        for(int i=0;i<shape[0];i++)
            for(int j=0;j<m;j++)
                B[i][j]=A[i*m+j];
        tensor<n,T> R(shape[0]);
        std::array<int,n-1> subshape;
        for(int i=1;i<n;i++)
            subshape[i-1]=shape[i];
        for(int i=0;i<shape[0];i++)
            R[i]=reshape<T,n-1>(B[i],subshape);
        return R;
    }
}


template<int n,bool is_inverse=false>
struct multidimensional_fft
{
    std::array<int,n> shape;
    multidimensional_fft(std::array<int,n> _shape):shape(std::move(_shape))
    {
    }
    using IC=std::complex<real>;
    tensor<n,IC> operator()(const tensor<n,IC>&T) const
    {
        tensor<n,IC> V(shape[0]);
        std::array<int,n-1> subshape;
        for(int i=1;i<n;i++)
            subshape[i-1]=shape[i];
        multidimensional_fft<n-1,is_inverse> subFFT(subshape);
        fast_fourier<is_inverse> FFT_1D(shape[0]);
        for(int i=0;i<shape[0];i++)
            V[i]=subFFT(T[i]);
        std::array<int,n-1> S;
        for(auto &s:S)
            s=0;
        std::vector<std::vector<IC>> R(shape[0]);
        do {
            std::vector<IC> Z;
            for(int i=0;i<shape[0];i++)
                Z.push_back(get<IC,n-1>(V[i],S));
            auto W=FFT_1D(Z);
            for(int i=0;i<shape[0];i++)
                R[i].push_back(W[i]);
            int k;
            for(k=0;k<n-1 && S[k] == subshape[k]-1;k++)
                S[k]=0;
            if(k<n-1)
                S[k]++;
        }while(std::any_of(S.begin(),S.end(),[](auto x)->bool{return x>0;}));
        tensor<n,IC> Y(shape[0]);
        for(int i=0;i<shape[0];i++)
            Y[i]=reshape<IC,n-1>(R[i],subshape);
        return Y;
    }
};

template<bool is_inverse>
struct multidimensional_fft<0,is_inverse>
{
    using IC=std::complex<real>;
    multidimensional_fft(std::array<int,0> shape){}
    tensor<0,IC> operator()(const tensor<0,IC> &O) const{
        return O;
    }
};

std::vector<IC> fast_multiply(std::vector<IC> x,std::vector<IC> y)
{
    int n=x.size(),m=y.size();
    int r=n+m;
    x.resize(r);
    y.resize(r);
    fast_fourier FFT(r);
    fast_fourier<true> IFFT(r);
    auto u=FFT(x),v=FFT(y);
    std::vector<IC> w(r);
    for(int i=0;i<r;i++)
        w[i]=u[i]*v[i];
    auto z=IFFT(w);
    for(auto &s:z)
        s/=r;
    z.resize(n+m-1);
    return z;
}

#endif //ACPC_PREPARATION_FFT_H

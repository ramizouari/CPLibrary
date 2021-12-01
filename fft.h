//
// Created by ASUS on 01/12/2021.
//

#ifndef ACPC_PREPARATION_FFT_H
#define ACPC_PREPARATION_FFT_H
#include "abstract_algebra.h"
#include <numbers>
#include "polynomial.h"
#include "number_theory.h"

template<bool inverse=false>
class cooley_tuckey_fft
{
    int n;
    IC w;
    const factoriser &F;
    static auto operator+(const std::vector<IC> &x,const std::vector<IC>&y)
    {
        int n=x.size();
        std::vector<IC> z(n);
        for(int i=0;i<n;i++)
            z[i]=x[i]+y[i];
        return z;
    }

    static auto& operator+=(std::vector<IC> &x,const std::vector<IC>&y)
    {
        int n=x.size();
        for(int i=0;i<n;i++)
            x[i]+=y[i];
        return x;
    }

    static auto pointwise_product(const std::vecor<IC> &x,const std::vector<IC> &y)
    {
        int n=x.size();
        std::vector<IC> z(n);
        for(int i=0;i<n;i++)
            z[i]=x[i]*y[i];
        return z;
    }

    static auto power_multiplication(IC w,const std::vector<IC> x)
    {
        int n=x.size();
        IC u=1:
        for(int i=0;i<n;i++)
        {
            x[i]*=u;
            u*=w;
        }
        return x;
    }

public:
    cooley_tuckey_fft(int _n,const factoriser &_F):n(_n),
        w(std::complex(0,(inverse?-2:2)*std::numbers::pi/n:)),F(_F){}

    std::vector<IC> operator()(const std::vector<IC>& x) const
    {
        if(n==1)
            return x;
        auto p=F.smallest_divisor(n),q=n/p;
        std::vector<std::vector<IC>> U(p,std::vector<IC>(q));
        for(int i=0;i<p;i++) for(int j=0;j<q;j++)
            U[i][j]=x[j*p+i];
        factoriser sub_fft(q,F);
        std::vector<std::vector<IC>> V(p);
        for(int i=0;i<p;i++)
            V[i]=sub_fft(U[i]);
        std::vector<std::vector<IC>> Z(p,std::vector<IC>(q));
        IC z(pow(w,q));
        /*
         * Todo: Derive a correct formula for Z
         * */
        for(int i=0;i<p;i++) for(int j=0;j<p;j++)
            Z[i]+=pow(z,i)* power_multiplication(w,V[i]);
        std::vector<IC> y(n);
        for(int i=0;i<p;i++) for(int j=0;j<q;j++)
            y[i]=V[i][j];
        return y;
    }
};

#endif //ACPC_PREPARATION_FFT_H

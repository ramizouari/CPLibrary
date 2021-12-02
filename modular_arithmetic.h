//
// Created by ASUS on 30/11/2021.
//
#ifndef __MODULAR__ARITHMETIC__
#define __MODULAR__ARITHMETIC__
#include <cstdint>
#include <utility>
#include "abstract_algebra.h"
#include <random>
#include <unordered_map>


template<integer m>
class cyclic
{
    integer n;
public:
    cyclic(int o=0):n((o+m)%m){}
    bool operator==(int O) const
    {
        return n==(m+O)%m;
    }

    bool operator!=(int O) const
    {
        return n!=(m+O)%m;
    }

    bool operator==(cyclic O) const
    {
        return n==O.n;
    }

    bool operator!=(cyclic O) const
    {
        return n!=O.n;
    }

    cyclic operator-() const
    {
        return cyclic(-n);
    }

    auto& operator+=(const cyclic &O)
    {
        n=(n+O.n)%m;
        return *this;
    }
    auto& operator-=(const cyclic &O)
    {
        n=(n+m-O.n)%m;
        return *this;
    }

    auto& operator*=(const cyclic &O)
    {
        n=(n*O.n)%m;
        return *this;
    }

    auto& operator/=(const cyclic &O)
    {
        return (*this)*=O.inv();
    }

    auto operator*(const cyclic &O) const
    {
        auto w=*this;
        return w*=O;
    }

    auto operator+(const cyclic &O) const
    {
        auto w=*this;
        return w+=O;
    }

    auto operator-(const cyclic &O) const
    {
        auto w=*this;
        return w-=O;
    }

    auto operator/(const cyclic &O) const
    {
        return (*this)*O.inv();
    }

    auto inv() const
    {
        return pow(*this,m-2);
    }
    auto pinv() const
    {
        return egcd(n,m).a;
    }

    auto& operator++()
    {
        return *this+=1;
    }

    auto& operator--()
    {
        return *this-=1;
    }

    explicit operator integer&()
    {
        return n;
    }

    explicit operator const integer&() const
    {
        return n;
    }
};

class d_cyclic
{
    integer n;
public:
    inline static integer m=1;
    d_cyclic(int o=0):n((o+m)%m){}
    bool operator==(int O) const
    {
        return n==(m+O)%m;
    }

    bool operator!=(int O) const
    {
        return n!=(m+O)%m;
    }

    bool operator==(d_cyclic O) const
    {
        return n==O.n;
    }

    bool operator!=(d_cyclic O) const
    {
        return n!=O.n;
    }

    auto& operator+=(const d_cyclic &O)
    {
        n=(n+O.n)%m;
        return *this;
    }
    auto& operator-=(const d_cyclic &O)
    {
        n=(n+m-O.n)%m;
        return *this;
    }

    auto& operator*=(const d_cyclic &O)
    {
        n=(n*O.n)%m;
        return *this;
    }

    auto inv() const
    {
        return pow(*this,m-2);
    }

    auto& operator/=(const d_cyclic &O)
    {
        return (*this)*=O.inv();
    }

    auto operator*(const d_cyclic &O) const
    {
        auto w=*this;
        return w*=O;
    }

    auto operator+(const d_cyclic &O) const
    {
        auto w=*this;
        return w+=O;
    }

    d_cyclic operator-() const
    {
        return d_cyclic(m-n);
    }

    auto operator-(const d_cyclic &O) const
    {
        auto w=*this;
        return w-=O;
    }

    auto operator/(const d_cyclic &O) const
    {
        return (*this)*O.inv();
    }

    auto pinv() const
    {
        return egcd(n,m).a;
    }

    auto& operator++()
    {
        return *this+=1;
    }

    auto& operator--()
    {
        return *this-=1;
    }

    explicit operator integer&()
    {
        return n;
    }

    explicit operator const integer&() const
    {
        return n;
    }
};

template<>
struct std::hash<d_cyclic>
{
    inline static std::random_device dev=std::random_device();
    inline static std::mt19937 g=std::mt19937(dev());
    inline static constexpr integer M=1e9+7;
    std::uniform_int_distribution<integer> d=std::uniform_int_distribution<integer>(1,M);
    integer a=d(g),b=d(g);
public:
    auto operator()(const d_cyclic &x) const
    {
        return (a*static_cast<integer>(x)+b)%M;
    }
};

template<int m>
struct std::hash<cyclic<m>>
{
    inline static std::random_device dev=std::random_device();
    inline static std::mt19937 g=std::mt19937(dev());
    inline static constexpr integer M=1e9+7;
    std::uniform_int_distribution<integer> d=std::uniform_int_distribution<integer>(1,M);
    integer a=d(g),b=d(g);
public:
    auto operator()(const cyclic<m> &x) const
    {
        return (a*static_cast<integer>(x)+b)%M;
    }
};

integer discrete_log(d_cyclic a,d_cyclic r)
{
    integer s=std::ceil(std::sqrt(d_cyclic::m));
    d_cyclic u=pow(a,s),w=1;
    std::unordered_map<d_cyclic,integer> mapper;
    for(integer i=0;i<=s;i++,w*=a)
        mapper[r*w]=i;
    w=u;
    for(integer i=1;i<=s;i++,w*=u)
        if(mapper.count(w))
            return i*s-mapper[w];
    return -1;
}

template<integer m>
integer discrete_log(cyclic<m> a,cyclic<m> r)
{
    static integer s=std::ceil(std::sqrt(m));
    cyclic<m> u=pow(a,s),w=1;
    std::unordered_map<d_cyclic,integer> mapper;
    for(integer i=0;i<=s;i++,w*=a)
        mapper[r*w]=i;
    w=u;
    for(integer i=1;i<=s;i++,w*=u)
        if(mapper.count(w))
            return i*s-mapper[w];
    return -1;
}

std::vector<integer> inverse_table(int n,int prime)
{
    std::vector<integer> I(n + 1);
    I[0] = I[1] = 1;
    for (int i = 2; i <= n; i++)
        I[i] = I[prime % i] *
                (prime - prime / i) % prime;
    return I;
}
#endif
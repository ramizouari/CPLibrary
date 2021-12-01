//
// Created by ASUS on 30/11/2021.
//
#include <cstdint>
#include <utility>
#include "abstract_algebra.h"


template<integer m>
class cyclic
{
    integer n;
public:
    cyclic(int o=0):n((o+m)%m){}
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
};

class d_cyclic
{
    integer n;
public:
    inline static integer m=1;
    d_cyclic(int o=0):n((o+m)%m){}
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
};
#include <iostream>
int main()
{
    auto [a,b,d]=egcd(17,25);
    std::cout << a << ' ' << b << ' ' << d;
}
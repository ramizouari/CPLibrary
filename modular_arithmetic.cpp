//
// Created by ASUS on 30/11/2021.
//
#include <cstdint>
#include <utility>

using integer=std::int64_t;
struct egcd_t
{
    integer a,b,d;
};

egcd_t egcd(integer a,integer b)
{
    if(a<b)
    {
        auto e = egcd(b, a);
        std::swap(e.a,e.b);
        return e;
    }
    integer q,s=1,t=0,tmp;
    q=a/b;
    tmp=b;
    b=a-b*q;
    a=tmp;
    while(b!=0)
    {
        s=s-b*q;
        t=t-b*q;
        q=a/b;
        tmp=b;
        b=a-b*q;
        a=tmp;
    }
    return {s,t,a};
}

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
};
#include <iostream>
int main()
{
    auto [a,b,d]=egcd(15,5);
    std::cout << a << ' ' << b << ' ' << d;
}
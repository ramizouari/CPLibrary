//
// Created by ramizouari on 01/12/2021.
//

#ifndef ACPC_PREPARATION_ABSTRACT_ALGEBRA_H
#define ACPC_PREPARATION_ABSTRACT_ALGEBRA_H

using natural = std::uint64_t;
using integer = std::int64_t;
using real = long double;

template<typename R>
R commutator(R a,R b)
{
    return a*b-b*a;
}

template<typename R>
R pow(R a, long long n)
{
    if(n==0)
        return 1;
    else if(n==1)
        return a;
    auto s=pow(a,n/2);
    return n%2?s*s*a:s*s;
}

template<typename R=integer>
struct egcd_t
{
    R a,b,d;
};

template<typename R>
egcd_t<R> egcd(R a,R b)
{
    if(a<b)
    {
        auto e = egcd(b, a);
        std::swap(e.a,e.b);
        return e;
    }
    integer q,s1=1,s2=0,t1=0,t2=1,tmp;
    while(b!=0)
    {
        q=a/b;
        tmp=s2;
        s2=s1-q*s2;
        s1=tmp;
        tmp=t2;
        t2=t1-q*t2;
        t1=tmp;
        tmp=b;
        b=a-b*q;
        a=tmp;
    }
    return {s1,t1,a};
}

#endif //ACPC_PREPARATION_ABSTRACT_ALGEBRA_H

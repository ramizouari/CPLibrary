//
// Created by ramizouari on 01/12/2021.
//

#ifndef ACPC_PREPARATION_ABSTRACT_ALGEBRA_H
#define ACPC_PREPARATION_ABSTRACT_ALGEBRA_H
#include <complex>
#include <functional>
using natural = std::uint64_t;
using integer = std::int64_t;
using real = long double;
using IR=real;
using IC= std::complex<IR>;

template<typename R>
R commutator(R a,R b)
{
    return a*b-b*a;
}


template<typename M,typename Endomorphism,typename G=typename M::base_field>
void apply_pointwise(const Endomorphism& f, M& u)
{
    if constexpr (std::is_same_v<G, M>)
        u = f(u);
    else for (auto& w : u)
        apply_pointwise<std::remove_reference_t<decltype(w)>, Endomorphism,G>(f, w);
}

template<typename M, typename BinaryOperation, typename G = typename M::base_field>
void apply_pointwise(const BinaryOperation& f, M& u, M a,M b)
{
    if constexpr (std::is_same_v<G, M>)
        u = f(a,b);
    else for (auto [p,s,t] : zip<M,M,M>(u,a,b))
        apply_pointwise<std::remove_reference_t<decltype(p)>, BinaryOperation, G>(f, p,s,t);
}

template<typename M, typename Endomorphism, typename G = typename M::base_field>
M pointwise_function(const Endomorphism& f, M u)
{
    apply_pointwise<M, Endomorphism,G>(f, u);
    return u;
}

template<typename M, typename BinaryOperation, typename G = typename M::base_field>
M pointwise_function(const BinaryOperation& f, M u,M a,M b)
{
    apply_pointwise<M, BinaryOperation, G>(f, u,a,b);
    return u;
}

template<typename R,typename ...StructureMetaData>
R pow(R a, long long n,StructureMetaData ... meta_info)
{
    if(n==0)
        return R(1,meta_info...);
    else if(n==1)
        return a;
    auto s=pow(a,n/2);
    return n%2?s*s*a:s*s;
}

template<typename R>
R gcd(R a,R b)
{
    if(a<b)
        std::swap(a,b);
    R q,tmp;
    while(b!=0)
    {
        q=a/b;
        tmp=b;
        b=a-b*q;
        a=tmp;
    }
    return a;
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
    R q,s1=1,s2=0,t1=0,t2=1,tmp;
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

template<typename R>
std::pair<R,R> bezout(R a, R b)
{
    auto [u,v,_]=egcd(a,b);
    return {u,v};
}

#endif //ACPC_PREPARATION_ABSTRACT_ALGEBRA_H

//
// Created by ramizouari on 01/12/2021.
//
#ifndef __NT_H__
#define __NT_H__
#include <vector>
#include "modular_arithmetic.h"

template<integer m>
struct uniform_hash
{
    using R=cyclic<m>;
    inline static std::random_device dev;
    std::mt19937_64 g;
    std::uniform_int_distribution<long long> d;
    R a,b;
    uniform_hash():g(dev()),d(1,m-1),a(d(g)),b(d(g)){}
    integer operator()(const integer &_x) const noexcept
    {
        R x(_x);
        return static_cast<integer>(a*x+b);
    }
};

template<integer n>
struct almost_uniform_hash
{
    using R=cyclic<n>;
    inline static std::random_device dev;
    std::mt19937_64 g;
    std::uniform_int_distribution<long long> d;
    R x;
    almost_uniform_hash():g(dev()),d(1,n-1),x(d(g)){}
    R operator()(const std::string &A) const noexcept
    {
        R r=0,w=1;
        for(auto a:A)
        {
            r=r+R(a)*w;
            w*=x;
        }
        return r;
    }
    std::vector<R> prefix_table(const std::string &A)const noexcept
    {
        std::vector<R> P;
        P.reserve(A.size()+1);
        P.emplace_back(0);
        R r=0,w=1;
        for(auto a:A)
        {
            r=r+R(a)*w;
            P.push_back(r);
            w*=x;
        }
        return P;
    }
};


template<long long n0,long long ...n>
class rabin_karp
{
    almost_uniform_hash<n0> H;
    rabin_karp<n...> RK;
    using IF=cyclic<n0>;
public:
    std::vector<int> potential_matches(const std::string &P,const std::string &T)
    {
        auto recursive_matches=RK.potential_matches(P,T);
        if(recursive_matches.empty())
            return {};
        std::vector<int> matches;
        int m=T.size(),r=P.size();
        auto A=H.prefix_table(T);
        auto R=H(P);
        std::vector<IF> X(m);
        X[0]=1;
        for(int i=1;i<m;i++)
            X[i]=X[i-1]*H.x;

        for(auto i:recursive_matches)
            if((A[i+r]-A[i])==(R*X[i]))
                matches.push_back(i);
        return matches;
    }

    int match(const std::string &P,const std::string &T)
    {
        auto matches=potential_matches(P,T);
        return matches.empty()?-1:matches.front();
    }

};

template<long long n>
class rabin_karp<n>
{
    almost_uniform_hash<n> H;
    using IF=cyclic<n>;
public:
    int match(const std::string &P,const std::string &T)
    {
        auto matches=potential_matches(P,T);
        return matches.empty()?-1:matches.front();
    }
    std::vector<int> potential_matches(const std::string &P,const std::string &T)
    {
        std::vector<int> matches;
        int m=T.size(),r=P.size();
        auto A=H.prefix_table(T);
        auto R=H(P);
        std::vector<IF> X(m);
        X[0]=1;
        for(int i=1;i<m;i++)
            X[i]=X[i-1]*H.x;
        for(int i=0;i<m-r+1;i++)
            if((A[i+r]-A[i])==(R*X[i]))
                matches.push_back(i);
        return matches;
    }
};

#endif
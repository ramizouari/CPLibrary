//
// Created by ramizouari on 01/12/2021.
//
#include "abstract_algebra.h"
#include "polynomial.h"

template<typename R>
class rational_extension
{
    R p,q;
    void reduce()
    {
        auto d=gcd(p,q);
        p/=d;
        q/=d;
    }
public:
    rational_extension(R _p=0,R _q=1):p(_p),q(_q)
    {
        reduce();
    }
    bool operator==(R a) const
    {
        return p==a*q;
    }
    bool operator!=(R a) const
    {
        return p!=a*q;
    }
    auto& operator+=(const rational_extension &o)
    {
        p=p*o.q+o.p*q;
        q*=o.q;
        reduce();
    }

    auto& operator-=(const rational_extension &o)
    {
        p=p*o.q-o.p*q;
        q*=o.q;
        reduce();
    }

    auto& operator*=(const rational_extension &o)
    {
        p*=o.p;
        q*=o.q;
        reduce();
    }

    auto operator+(const rational_extension &o) const
    {
        auto r=*this;
        return r+=o;
    }

    auto operator-(const rational_extension &o) const
    {
        auto r=*this;
        return r-=o;
    }

    auto operator*(const rational_extension &o) const
    {
        auto r=*this;
        return r*=o;
    }

    auto operator-() const
    {
        return rational_extension(-p,q);
    }
};

template<typename R,int nilpotence>
class nilpotent_extension
{
    std::vector<R> p;
    void reduce()
    {
        while(!p.empty() && p.back()==0)
            p.pop_back();
    }
public:
    nilpotent_extension(R k=0):p(1,k)
    {
        reduce();
    }
    nilpotent_extension(const std::vector<R> &_p)
    {
        int n=_p.size();
        p.resize(std::min(n,nilpotence));
        for(int i=0;i<std::min(n,nilpotence);i++)
            p[i]=_p[i];
        reduce();
    }
    auto& operator+=(const nilpotent_extension &O)
    {
        p.resize(std::max(p.size()),O.p.size());
        for(int i=0;i<O.size();i++)
            p[i]+=O.p[i];
        reduce();
        return *this;
    }

    auto& operator-=(const nilpotent_extension &O)
    {
        p.resize(std::max(p.size()),O.p.size());
        for(int i=0;i<O.size();i++)
            p[i]-=O.p[i];
        reduce();
        return *this;
    }

    auto operator*(const nilpotent_extension &O) const
    {
        nilpotent_extension q;
        int n=p.size()-1,m=O.p.size()-1;
        q.p.resize(std::min(n+m+1,nilpotence));
        for(int i=0;i<n;i++) for(int j=0;j<m && (i+j) < nilpotence;j++)
            q.p[i+j]+=p[i]*O.p[j];
        reduce();
        return q;
    }

    auto& operator*=(const nilpotent_extension &O)
    {
        auto q=(*this)*O;
        p.swap(q.p);
        return *this;
    }

    auto operator+(const nilpotent_extension &O) const
    {
        auto q=*this;
        return q+=O;
    }

    auto operator-(const nilpotent_extension &O) const
    {
        auto q=*this;
        return q-=O;
    }
};

template<typename R,int idempotence>
class idempotent_extension
{
    std::vector<R> p;
    void reduce()
    {
        while(!p.empty() && p.back()==0)
            p.pop_back();
    }
public:
    idempotent_extension(R k=0):p(1,k)
    {
        reduce();
    }
    idempotent_extension(const std::vector<R> &_p)
    {
        int n=_p.size();
        p.resize(std::min(n,idempotence));
        for(int i=0;i<n;i++)
            p[i]=_p[std::min(i,idempotence-1)];
        reduce();
    }
    auto& operator+=(const idempotent_extension &O)
    {
        p.resize(std::max(p.size()),O.p.size());
        for(int i=0;i<O.size();i++)
            p[i]+=O.p[i];
        reduce();
        return *this;
    }

    auto& operator-=(const idempotent_extension &O)
    {
        p.resize(std::max(p.size()),O.p.size());
        for(int i=0;i<O.size();i++)
            p[i]-=O.p[i];
        reduce();
        return *this;
    }

    auto operator*(const idempotent_extension &O) const
    {
        idempotent_extension q;
        int n=p.size()-1,m=O.p.size()-1;
        q.p.resize(std::min(n+m+1,idempotence));
        for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                q.p[std::min(i+j,idempotence-1)]+=p[i]*O.p[j];
        reduce();
        return q;
    }

    auto& operator*=(const idempotent_extension &O)
    {
        auto q=(*this)*O;
        p.swap(q.p);
        return *this;
    }

    auto operator+(const idempotent_extension &O) const
    {
        auto q=*this;
        return q+=O;
    }

    auto operator-(const idempotent_extension &O) const
    {
        auto q=*this;
        return q-=O;
    }
};

template<typename R>
class ring_extension
{
    polynomial<R> p;
public:
    void reduce()
    {
        if(p.degree()>=q.degree())
            p=p%q;
        p.reduce();
    }

    inline static polynomial<R> q;
    ring_extension(R k=0):p(1,k)
    {
        reduce();
    }
    ring_extension(const std::vector<R> &_p):p(_p)
    {
        reduce();
    }
    auto& operator+=(const ring_extension &O)
    {
        p+=O.p;
        return *this;
    }

    auto& operator-=(const ring_extension &O)
    {
        p-=O.p();
        return *this;
    }

    auto operator*(const ring_extension &O) const
    {
        auto q=*this;
        return q*=O;
    }

    auto& operator*=(const ring_extension &O)
    {
        p*=O.p;
        reduce();
        return *this;
    }

    auto operator+(const ring_extension &O) const
    {
        auto q=*this;
        return q+=O;
    }

    auto operator-(const ring_extension &O) const
    {
        auto q=*this;
        return q-=O;
    }
};
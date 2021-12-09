//
// Created by ramizouari on 01/12/2021.
//
#ifndef __POLYNOMIAL__H__
#define __POLYNOMIAL__H__
#include <vector>
#include <map>

template<typename R>
class polynomial
{
    std::vector<R> p;

public:
    void reduce()
    {
        while(!p.empty() && p.back()==R(0))
            p.pop_back();
    }
    polynomial(R k=0):p(1,k)
    {
        reduce();
    }
    polynomial(std::vector<R> _p):p(std::move(_p))
    {
        reduce();
    }

    int degree() const
    {
        return p.size()-1;
    }

    auto& operator+=(const polynomial &O)
    {
        int r=std::max(p.size(),O.p.size());
        p.resize(r);
        for(int i=0;i<O.p.size();i++)
            p[i]+=O.p[i];
        reduce();
        return *this;
    }

    auto& operator-=(const polynomial &O)
    {
        int r=std::max(p.size(),O.p.size());
        p.resize(r);
        for(int i=0;i<O.p.size();i++)
            p[i]-=O.p[i];
        reduce();
        return *this;
    }

    polynomial operator*(const polynomial &O) const
    {
        if(O.p.empty() || p.empty())
            return polynomial(0);
        int n=degree(),m=O.degree();
        polynomial q;
        q.p.resize(n+m+1);
        for(int i=0;i<=n;i++) for(int j=0;j<=m;j++)
            q.p[i+j]+=p[i]*O.p[j];
        q.reduce();
        return q;
    }

    auto& operator*=(const polynomial &O)
    {
        auto r=(*this)*O;
        p.swap(r.p);
        return *this;
    }

    auto operator+(const polynomial &O) const
    {
        auto r=*this;
        return r+=O;
    }

    auto operator-(const polynomial &O) const
    {
        auto r=*this;
        return r-=O;
    }

    auto operator-() const
    {
        auto r=*this;
        for(auto &s:r.p)
            s=-s;
        return r;
    }

    auto operator*=(R a)
    {
        if(a==0)
            p.clear();
        else for(auto& s:p)
            s*=a;
        reduce();
        return *this;
    }

    bool operator==(R a) const
    {
        if(a==0)
            return p.empty();
        else return degree() == 0 && p.front() == a;
    }

    bool operator!=(R a) const
    {
        return !(*this==a);
    }

    auto& operator+=(R a)
    {
        return *this+=polynomial({a});
    }

    auto& operator-=(R a)
    {
        return *this-=polynomial({a});
    }

    auto operator+(R a) const
    {
        auto q=*this;
        return q+=a;
    }

    auto operator-(R a) const
    {
        auto q=*this;
        return q-=a;
    }

    bool operator<(const polynomial &O) const
    {
        return degree() < O.degree();
    }

    auto& operator/=(R k)
    {
        for(auto &s:p)
            s/=k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto q=*this;
        return q/=k;
    }

    std::pair<polynomial,polynomial> euclidean_division(const polynomial &O) const
    {
        if(degree() < O.degree())
            return {R(0),*this};
        polynomial q,r=*this;
        int n=degree(),m=O.degree();
        q.p.resize(n-m+1);
        for(int i=n;i>=m;i--)
        {
            q.p[i-m]=r[i]/O.p[m];
            for(int j=0;j<=m;j++)
                r.p[i+j-m]-=q.p[i-m]*O.p[j];
        }
        r.p.resize(m);
        return {q,r};
    }

    polynomial operator/(const polynomial &O) const
    {
        return euclidean_division(O).first;
    }

    polynomial operator%(const polynomial &O) const
    {
        return euclidean_division(O).second;
    }

    auto& operator/=(polynomial &O)
    {
        p.swap(((*this)/O).p);
        return *this;
    }

    auto &operator[](int k)
    {
        return p[k];
    }

    const auto& operator[](int k) const
    {
        return p[k];
    }

    template<typename H>
    H operator()(H a)
    {
        H r=0;
        for(int i=degree();i>=0;i--)
            r=r*a+p[i];
        return r;
    }

    auto begin()
    {
        return p.begin();
    }

    auto end()
    {
        return p.end();
    }

    explicit operator std::vector<R>&()
    {
        return p;
    }

    explicit operator const std::vector<R>&() const
    {
        return p;
    }
};

template<typename R>
polynomial<R> operator*(R a,const polynomial<R> &p)
{
    auto q=p;
    return q*=a;
}

template<typename R>
const polynomial<R> X=polynomial<R>(std::vector<R>{0,1});

template<typename R>
class sparse_polynomial
{
    std::map<int,R> p;
    void reduce()
    {
        std::vector<int> to_del;
        for(auto [k,x]:p)
            if(x==0)
                to_del.push_back(k);
        for(auto k:to_del)
            p.erase(k);
    }
public:
    sparse_polynomial(R k=0)
    {
        p[0]=k;
        reduce();
    }
    sparse_polynomial(const std::vector<R> &_p)
    {
        for(int i=0;i<_p.size();i++)
            p[i]=_p[i];
        reduce();
    }

    int degree() const
    {
        return p.empty()?-1:p.rbegin()->first;
    }

    auto& operator+=(sparse_polynomial O)
    {
        for(const auto& [k,s]:O.p)
        {
            p[k] += O.p[k];
            if(p[k]==0)
                p.erase(k);
        }
        return *this;
    }

    auto& operator-=(sparse_polynomial O)
    {
        for(const auto& [k,s]:O.p)
        {
            p[k] -= O.p[k];
            if(p[k]==R(0))
                p.erase(k);
        }
        return *this;
    }

    auto operator*(const sparse_polynomial &O) const
    {

        sparse_polynomial q;
        for(auto [i,u]:p) for(auto [j,v]:O.p)
        {
            q.p[i+j]+=u*v;
            if(q.p[i+j]==0)
                q.p.erase(i+j);
        }
        return q;
    }

    auto& operator*=(const sparse_polynomial &O)
    {
        auto r=(*this)*O;
        p.swap(r.p);
        return *this;
    }

    auto operator+(const sparse_polynomial &O) const
    {
        auto r=*this;
        return r+=O;
    }

    auto operator-(const sparse_polynomial &O) const
    {
        auto r=*this;
        return r-=O;
    }

    auto operator-() const
    {
        auto r=*this;
        for(auto &s:r.p)
            s=-s;
        return r;
    }

    auto operator*=(R a)
    {
        if(a==R(0))
            p.clear();
        else for(auto& s:p)
                s*=a;
        reduce();
        return *this;
    }

    auto& operator+=(R a)
    {
        return *this+=sparse_polynomial({a});
    }

    auto& operator-=(R a)
    {
        return *this+=sparse_polynomial({a});
    }

    auto operator+(R a) const
    {
        auto q=*this;
        return q+=a;
    }

    auto operator-(R a) const
    {
        auto q=*this;
        return q-=a;
    }

    auto& operator/=(R k)
    {
        for(auto &s:p)
            s/=k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto q=*this;
        return q/=k;
    }

    auto &operator[](int k)
    {
        return p[k];
    }

    const auto& operator[](int k) const
    {
        return p[k];
    }

    template<typename H>
    H operator()(H a)
    {
        H r=0,u=1,i=0;
        for(auto [k,x]:p)
        {
            u*=pow(a,k-i);
            r+=u*x;
            i=k;
        }
    }
};

template<typename R>
sparse_polynomial<R> Z=sparse_polynomial<R>({0,1});

template<typename R>
sparse_polynomial<R> operator*(R a,const sparse_polynomial<R> &p)
{
    auto q=p;
    return q*=a;
}

template<typename R>
polynomial<R> newton_interpolation(const std::vector<R> &x,const std::vector<R> &y)
{
    int n=x.size()-1;
    std::vector<std::vector<R>> d(n+1,std::vector<R>(n+1));
    for(int i=0;i<=n;i++)
        d[i][i]=y[i];
    for(int r=1;r<=n;r++) for(int i=0;i+r<=n;i++)
        d[i][i+r]=(d[i+1][i+r]-d[i][i+r-1])/(x[i+r]-x[i]);
    polynomial<R> p,u=R(1);
    for(int i=0;i<=n;i++) {
        p +=d[0][i]*u;
        u*=(X<R>-x[i]);
    }
    return p;
}

template<typename R>
polynomial<R> karatsuba_multiplication(const polynomial<R> &p,const polynomial<R> &q)
{
    constexpr int L=75;
    if(std::min(p.degree(),q.degree())<=L)
        return p*q;
    polynomial<R> a1,b1,a2,b2;
    int n=p.degree(),m=q.degree(),r=std::max(n,m)+1;
    std::vector<R> &u1=static_cast<std::vector<R>&>(a1),&u2=static_cast<std::vector<R>&>(a2),
    &v1=static_cast<std::vector<R>&>(b1),&v2=static_cast<std::vector<R>&>(b2);
    u1.resize(std::min(n+1,r/2));
    u2.resize(std::min(m+1,r/2));
    v1.resize(std::max(n+1-r/2,0));
    v2.resize(std::max(m+1-r/2,0));
    for(int i=0;i<u1.size();i++)
        u1[i]=p[i];
    for(int i=0;i<u2.size();i++)
        u2[i]=q[i];
    for(int i=0;i<v1.size();i++)
        v1[i]=p[i+r/2];
    for(int i=0;i<v2.size();i++)
        v2[i]=q[i+r/2];
    polynomial<R> r1= karatsuba_multiplication(a1,a2),
        r3= karatsuba_multiplication(b1,b2),
        t=karatsuba_multiplication(a1+b1,a2+b2),
        r2=t-r1-r3;
    polynomial<R> h;
    int s=r-r%2;
    auto &c=static_cast<std::vector<R>&>(h);
    c.resize(n+m+1);
    for(int i=0;i<=r1.degree();i++)
        c[i]+=r1[i];
    for(int i=0;i<=r2.degree();i++)
        c[i+r/2]+=r2[i];
    for(int i=0;i<=r3.degree();i++)
        c[i+s]+=r3[i];
    return h;
}
#endif
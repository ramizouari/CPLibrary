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
        while(!p.empty() && p.back()==0)
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
        return *this+=polynomial({a});
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
            r=(r+p[i])*a;
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
};

template<typename R>
polynomial<R> operator*(R a,const polynomial<R> &p)
{
    auto q=p;
    return q*=a;
}

template<typename R>
const polynomial<R> X=polynomial<R>({0,1});

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
        return *this+=polynomial({a});
    }

    auto& operator-=(R a)
    {
        return *this+=polynomial({a});
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
    polynomial<R> p,u=1;
    for(int i=0;i<=n;i++) {
        p +=d[0][i]*u;
        u*=(X<R>-x[i]);
    }
    return p;
}
#endif
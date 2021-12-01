//
// Created by ASUS on 01/12/2021.
//

#include <vector>

template<typename R>
class polynomial
{
    std::vector<R> p;
    void reduce()
    {
        while(!p.empty() && p.back()==0)
            p.pop_back();
    }
public:
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
    }

    auto& operator-=(const polynomial &O)
    {
        int r=std::max(p.size(),O.p.size());
        p.resize(r);
        for(int i=0;i<O.p.size();i++)
            p[i]-=O.p[i];
    }

    auto operator*(const polynomial &O) const
    {
        if(O.p.empty() || p.empty())
            return polynomial(0);
        int n=degree(),m=O.degree();
        polynomial q;
        q.p.resize(n+m);
        for(int i=0;i<=n;i++) for(int j=0;j<=m;j++)
            q[i+j]+=p[i]*p[j];
        q.reduce();
        return q;
    }

    auto& operator*=(const polynomial &O) const
    {
        auto r=*this;
        return r*=O;
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
        return degree()==0 && (p.front()==a);
    }

    bool operator!=(R a) const
    {
        return degree()==-1 || p.front()!=a;
    }

    auto& operator+=(R a)
    {
        return *this+=polynomial({a});
    }

    auto& operator-=(R a)
    {
        return *this+=polynomial({a});
    }

    auto operator+(R a)
    {
        auto q=*this;
        return q+=a;
    }

    auto operator-(R a)
    {
        auto q=*this;
        return q-=a;
    }
    auto operator/(const polynomial &O) const
    {
        if(degree() < O.degree())
            return 0;
        polynomial q,r=*this;
        int n=degree(),m=O.degree();
        q.p.resize(n-m+1);
        for(int i=n;i>=m;i--)
        {
            q.p[i-m]=p[i]/O.p[m];
            for(int j=i-m;j<=r.degree();j++)
                r.p[j+m]-=q.p[j]*O.p[m];
        }
        return q;
    }

    auto& operator/=(polynomial &O)
    {
        p.swap(((*this)/O).p);
        return *this;
    }

    template<typename H>
    H operator()(H a)
    {
        H r=0;
        for(int i=degree();i>=0;i--)
            r=r+a*p[i];
        return r;
    }
};

template<typename R>
polynomial<R> newton_interpolation(const std::vector<R> &x,const std::vector<R> &y)
{
    int n,m;

}
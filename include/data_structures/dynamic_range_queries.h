//
// Created by ramizouari on 25/10/23.
//

#ifndef CPLIBRARY_DYNAMIC_RANGE_QUERIES_H
#define CPLIBRARY_DYNAMIC_RANGE_QUERIES_H
#include "range_queries.h"
#include <memory>

template<typename T>
struct sparse_array<T,void>
{
    int n,h;
    std::vector<std::vector<T>> S;
    std::shared_ptr<binary_operation<T>> F;
public:
    sparse_array(const std::vector<T>&A, std::shared_ptr<binary_operation<T>> _F):n(bit_ceil(A.size())),h(bit_log(n)),S(h+1),F(_F)
    {
        int r=1;
        for(int i=h;i>=0;i--,r*=2)
            S[i].resize(n-r+1,F->neutral_element());
        for(int i=0;i<A.size();i++)
            S[h][i]=A[i];
        r=1;
        for(int i=h-1;i>=0;i--,r*=2) for(int j=0;j<=n-2*r;j++)
                S[i][j]=(*F)(S[i+1][j],S[i+1][j+r]);
    }

    T query(int l,int r) const
    {
        if(l>=r)
            return F->neutral_element();
        auto d=r-l;
        auto s=bit_floor(d);
        auto b=bit_log(s);
        return (*F)(S[h-b][l],S[h-b][r-s]);
    }
};

template<typename R>
struct segment_tree<R,void>
{
    std::vector<std::vector<R>> S;
    std::vector<R> A;
    int n,h;
    std::shared_ptr<binary_operation<R>> F;
    segment_tree(const std::vector<R> &_A, std::shared_ptr<binary_operation<R>> _F):A(_A),F(_F)
    {
        n=bit_ceil(A.size());
        A.resize(n,F->neutral_element());
        int m=n;
        h=0;
        while(m)
        {
            m/=2;
            h++;
        }
        S.resize(h);
        for(int i=0;i<h;i++)
            S[i].resize(1<<i);
        build();
    }

    void update(int i,R u)
    {
        A[i]=u;
        S[h-1][i]=u;
        int m=h-2;
        i/=2;
        while(m>=0)
        {
            S[m][i]=(*F)(S[m+1][2*i],S[m+1][2*i+1]);
            m--;
            i/=2;
        }
    }

    R query(int l,int r)
    {
        return query(std::max(l,0),std::min(r,n),0,n,0);
    }
private:
    void build()
    {
        for(int i=0;i<n;i++)
            S.back()[i]=A[i];
        for(int i=h-2;i>=0;i--) for(int k=0;k<(1<<i);k++)
                S[i][k]=(*F)(S[i+1][2*k],S[i+1][2*k+1]);
    }
    R query(int l,int r,int a,int b,int depth)
    {
        if(l>=r)
            return F->neutral_element();
        if(l==a && r==b)
            return S[depth][l>>(h-1-depth)];
        int mid=(a+b)/2;
        if(mid>r)
            return query(l,r,a,mid,depth+1);
        else if(mid<l)
            return query(l,r,mid,b,depth+1);
        else
            return (*F)(query(l,mid,a,mid,depth+1),query(mid,r,mid,b,depth+1));
    }
};

template<typename R>
struct prefix_array<R,void>
{
    std::vector<R> A;
    std::vector<R> P;
    std::shared_ptr<binary_operation<R>> F;
    std::shared_ptr<invertible_operation<R>> I;
    prefix_array(const std::vector<R> &_A, std::shared_ptr<binary_operation<R>> _F, std::shared_ptr<invertible_operation<R>> _I):A(_A),P(_A.size()+1),F(_F),I(_I)
    {
        P[0]=F->neutral_element();
        for(int i=0;i<A.size();i++)
            P[i+1]=(*F)(P[i],A[i]);
    }

    R query(int l,int r)
    {
        return (*F)(I->inv(P[l]),P[r]);
    }

    void update(int i,R u)
    {
        A[i]=u;
        for(int j=i+1;j<P.size();j++)
            P[j]=(*F)(P[j-1],A[j-1]);
    }
};


template<typename T>
struct fenwick_tree<T,void> {
    int n;
    std::shared_ptr<binary_operation<T>> F;
    std::shared_ptr<invertible_operation<T>> I;
    std::vector<T> bit;

    fenwick_tree(int _n,std::shared_ptr<binary_operation<T>> _F, std::shared_ptr<invertible_operation<T>> _I):n(_n),F(_F),I(_I),bit(n,F->neutral_element()){}
    fenwick_tree(const std::vector<T> &X) : fenwick_tree(X.size())
    {
        for(int i=0;i<n;i++)
            update(i,X[i]);
    }
    T sum(int x) {
        if(x<0)
            return F->neutral_element();
        T ret = F->neutral_element();
        for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
            ret = (*F)(ret,bit[i]);
        return ret;
    }

    T query(int a,int b)
    {
        return F(F.inv(sum(a-1)),sum(b));
    }

    T sum(int a,int b)
    {
        return query(a,b);
    }

    void add(int x, T delta) {
        for (int i = x; i < n; i = i | (i + 1))
            bit[i] = (*F)(bit[i], delta);
    }

    void update(int x, T delta) {
        add(x,(*F)(I->inv(sum(x,x)),delta));
    }
};


#endif //CPLIBRARY_DYNAMIC_RANGE_QUERIES_H

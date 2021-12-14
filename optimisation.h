//
// Created by ASUS on 01/12/2021.
//

#ifndef __OPTIMISATION__
#define __OPTIMISATION__
#include "linear_algebra.h"
#include "abstract_algebra.h"
#include "topology.h"
#include <functional>

template<int n>
s_vector<real,n> point_wise_divide(const s_vector<real,n> &x,const s_vector<real,n> &y)
{
    s_vector<real,n> z;
    for(int i=0;i<n;i++)
        z[i]=x[i]/y[i];
    return z;
}

template<int n,int m>
s_vector<real,n> get_column(const s_matrix<real,n,m> &A,int k)
{
    s_vector<real,n> C;
    for(int i=0;i<n;i++)
        C[i]=A[i][k];
    return C;
}

template<int n>
real max(const s_vector<real,n> &x)
{
    real r=x[0];
    for(auto s:x)
        r=std::max(r,s);
    return r;
}

template<int n>
int argmax(const s_vector<real,n> &x)
{
    int k=0;
    for(int i=1;i<n;i++) if(x[i]>x[k])
            k=i;
    return k;
}

template<int n, int m>
s_vector<real,m> simplex(
        const s_vector<real,m>& _Z,
        const s_matrix<real,n, m>& _A,
        s_vector<real,n> b)
{
    s_vector<real, n + m > Z;
    for (int i = 0; i < m; i++)
        Z[i] = _Z[i];
    s_matrix<real,n, n + m>A;
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
            A[i][j] = _A[i][j];
    for (int i = m; i < n + m; i++)
        A[i - m][i] = 1;
    while (max(Z) > 0)
    {
        auto q = argmax(Z);
        int p = -1;
        auto c = point_wise_divide(b,get_column(A,q));
        for (int k = 0; k < n; k++)
            if (A[k][q] > 0 && c[k] >= 0 && (p == -1 || c[k] < c[p]))
                p = k;
        if (p == -1)
            break;
        for (int i = 0; i < n; i++) if (i != p)
            {
                auto k=A[i][q];
                b[i] -= (k / A[p][q]) * b[p];
                for(int j=0;j<n+m;j++)
                    A[i][j] -= (k / A[p][q]) * A[p][j];
                A[i][q] = 0;

            }
        auto k=Z[q];
        for(int j=0;j<n+m;j++)
            Z[j] -= (k / A[p][q]) * A[p][j];
        Z[q] = 0;
    }
    s_vector<real,n + m> h;
    for (int i = 0; i < n; i++)
        h[i] = b[i];
    s_matrix<real,n + m, n + m> Q;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n + m; j++)
            Q[i][j] = A[i][j];
    int r = 0;
    for (int i = 0; i < n + m ; i++)
        if (Z[i] < 0)
        {
            Q[n + r++][i] = 1;
            for (int j = 0; j < n; j++)
                A[j][i] = 0;
        }
    s_vector<real,m> d;
    auto g = Q.solve(h);
    for (int i = 0; i < m; i++)
        d[i] = g[i];
    return d;
}


d_vector<real> point_wise_divide(const d_vector<real> &x,const d_vector<real> &y)
{
    int n=x.dim();
    d_vector<real> z(v_shape{n});
    for(int i=0;i<n;i++)
        z[i]=x[i]/y[i];
    return z;
}

d_vector<real> get_column(const d_matrix<real> &A,int k)
{
    int n=A.row_dim(),m=A.col_dim();
    d_vector<real> C(v_shape{m});
    for(int i=0;i<n;i++)
        C[i]=A[i][k];
    return C;
}

real max(const d_vector<real> &x)
{
    real r=x[0];
    for(auto s:x)
        r=std::max(r,s);
    return r;
}

int argmax(const d_vector<real> &x)
{
    int k=0,n=x.dim();
    for(int i=1;i<n;i++) if(x[i]>x[k])
        k=i;
    return k;
}

d_vector<real> simplex(
        const d_vector<real>& _Z,
        const d_matrix<real>& _A,
        d_vector<real> b)
{
    int n=b.dim(),m=_Z.dim();
    d_vector<real> Z(v_shape{n+m});
    for (int i = 0; i < m; i++)
        Z[i] = _Z[i];
    d_matrix<real>A(0,m_shape{n,n+m});
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
            A[i][j] = _A[i][j];
    for (int i = m; i < n + m; i++)
        A[i - m][i] = 1;
    while (max(Z) > 0)
    {
        auto q = argmax(Z);
        int p = -1;
        auto c = point_wise_divide(b,get_column(A,q));
        for (int k = 0; k < n; k++)
            if (A[k][q] > 0 && c[k] >= 0 && (p == -1 || c[k] < c[p]))
                p = k;
        if (p == -1)
            break;
        for (int i = 0; i < n; i++) if (i != p)
            {
                auto k=A[i][q];
                b[i] -= (k / A[p][q]) * b[p];
                for(int j=0;j<n+m;j++)
                    A[i][j] -= (k / A[p][q]) * A[p][j];
                A[i][q] = 0;

            }
        auto k=Z[q];
        for(int j=0;j<n+m;j++)
            Z[j] -= (k / A[p][q]) * A[p][j];
        Z[q] = 0;
    }
    d_vector<real> h(v_shape{n+m});
    for (int i = 0; i < n; i++)
        h[i] = b[i];
    d_matrix<real> Q(0,m_shape{n+m,n+m});
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n + m; j++)
            Q[i][j] = A[i][j];
    int r = 0;
    for (int i = 0; i < n + m ; i++)
        if (Z[i] < 0)
        {
            Q[n + r++][i] = 1;
            for (int j = 0; j < n; j++)
                A[j][i] = 0;
        }
    d_vector<real> d(v_shape{m});
    auto g = Q.solve(h);
    for (int i = 0; i < m; i++)
        d[i] = g[i];
    return d;
}

template<typename E,typename F,typename R>
class derivator
{
    real eps;
public:
    derivator(real _eps=1e-7):eps(_eps){}
    template<typename ... StructureMetaData>
    R jacobian(const std::function<F(E)>&f,E a,StructureMetaData ...meta_info) const
    {
        R J(0,meta_info...);
        for(int i=0;i<J.col_dim();i++)
        {
            a[i]-=eps;
            auto z2=f(a);
            a[i]+=2*eps;
            auto z1=f(a);
            auto du=(z2-z1)/(2*eps);
            for(int j=0;j<J.row_dim();j++)
                J[i][j]=du[j];
            a[i]-=eps;
        }
    }
};

template<typename E>
class derivator<E,real,E>
{
    real eps;
public:
    derivator(real _eps=1e-7):eps(_eps){}

    E gradient(const std::function<real(E)>&f,E a) const
    {
        E grad(a);
        for(int i=0;i<grad.dim();i++)
        {
            a[i]+=eps;
            auto z2=f(a);
            a[i]-=2*eps;
            auto z1=f(a);
            a[i]+=eps;
            grad[i]=(z2-z1)/(2*eps);
        }
        return grad;
    }
};

template<>
class derivator<real,real,real>
{
    real eps;
public:
    derivator(real _eps=1e-7):eps(_eps){}
    real derivative(const std::function<real(real)>&f,real a) const
    {
        return (f(a+eps)-f(a-eps))/(2*eps);
    }

    real gradient(const std::function<real(real)>&f,real a) const
    {
        return derivative(f,a);
    }
};

template<typename E,typename Norm=L2_inner_product<E>>
class gradient_descent
{
    inline static constexpr Norm N=Norm();
protected:
    E x0;
    real p=.1;
    real eps=1e-8;
    derivator<E,real,E>& D;
public:
    gradient_descent(const E& _x0,derivator<E,real,E> &d):x0(_x0),D(d) {}
    void set_seed(const E& _x0) { x0 = _x0; }
    E argmin(const std::function<real(E)>& f) const
    {
        E x = x0;
        for (; N.norm(D.gradient(f, x)) > eps; x -= p * D.gradient(f, x));
        return x;
    }
};

template<typename E,typename InnerProduct=L2_inner_product<E>>
class barzilai_borwein_gradient_descent
{
    E s;
    E x0;
    real p=.1;
    real eps=1e-8;
    derivator<E,real,E>& D;
    inline static constexpr InnerProduct B = InnerProduct();
public:
    barzilai_borwein_gradient_descent(const E& _x0, derivator<E, real,E>& d, real _p):x0(_x0),D(d),p(_p){}

    E argmin(const std::function<real(E)>& f)
    {
        this->p = 0.1;
        s = this->x0;
        E x = s- this->p*this->D.gradient(f, s);
        for (; B.norm(this->D.gradient(f, x)) > this->eps; x -= this->p * this->D.gradient(f, x))
        {
            update_rate(f, x);
            s = x;
        }
        return x;
    }

    virtual void update_rate(const std::function<real(E)>& f, const E& x)
    {
        auto L = this->D.gradient(f, x) - this->D.gradient(f, s);
        this->p = B.inner_product(L,x - s) / B.inner_product(L,L);
    }
};

#endif

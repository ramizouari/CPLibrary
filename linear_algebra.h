//
// Created by ramizouari on 30/11/2021.
//
#ifndef __LINEAR__ALGEBRA__
#define __LINEAR__ALGEBRA__
#include <vector>
#include <array>
#include "polynomial.h"

struct v_shape
{
    int n;
};

/*
 * Dynamic Vector:
 * It is the union of R^k for all k
 * - Addition between 2 vectors are defined with respect to the first vector's shape
 * - for all k, the set of all vectors of shape k is a vector space
 * @Requirements
 * R is a commutative ring
 * */
template<typename R>
class d_vector
{
    std::vector<R> u;
public:
    using base_field=R;
    using base_ring=R;
    inline static int n=0;
    d_vector():u(n){}
    d_vector(std::vector<R> _u):u(std::move(_u)){}
    d_vector(v_shape shape):u(shape.n){}
    auto dim() const
    {
        return u.size();
    }

    auto& operator[](int k)
    {
        return u[k];
    }

    const auto& operator[](int k) const
    {
        return u[k];
    }

    auto& operator+=(const d_vector &o)
    {
        for(int i=0;i<dim();i++)
            u[i]+=o.u[i];
        return *this;
    }

    auto& operator-=(const d_vector &o)
    {
        for(int i=0;i<dim();i++)
            u[i]-=o.u[i];
        return *this;
    }

    auto& operator*=(R k)
    {
        for(auto &s:u)
            s*=k;
        return *this;
    }

    auto operator+(const d_vector &o) const
    {
        auto v=*this;
        return v+=o;
    }

    auto operator-(const d_vector &o) const
    {
        auto v=*this;
        return v-=o;
    }

    auto operator-() const
    {
        auto v=*this;
        for(auto &s:v.u)
            s=-s;
        return v;
    }

    auto& operator/=(R k)
    {
        for(auto &s:u)
            s/=k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto v=*this;
        return v/=k;
    }

    auto begin() {
        return u.begin();
    }

    auto begin() const
    {
        return u.cbegin();
    }

    auto end()
    {
        return u.end();
    }

    auto end() const
    {
        return u.cend();
    }
};

template<typename R>
auto operator*(const R&k,const d_vector<R>& u)
{
    auto v=u;
    return v*=k;
}

/*
 * Static Vector:
 * It is a member of an R-vector space E where dim(E)= n
 * @Requirements
 * R is a commutative ring. Formally, E is a R-module, and it is a vector space only if R is a field; the name s_vector is used
 * for consistency with the computer science's name.
 */

template<typename R,int n>
class s_vector
{
    std::array<R,n> u;
public:
    using base_field=R;
    using base_ring=R;
    inline static constexpr int dim()
    {
        return n;
    }

    s_vector()
    {
        for(int i=0;i<n;i++)
            u[i]=0;
    }

    s_vector(std::array<R,n>_u):u(std::move(_u)){}

    auto& operator[](int k)
    {
        return u[k];
    }

    const auto& operator[](int k) const
    {
        return u[k];
    }

    auto& operator+=(const s_vector &o)
    {
        auto r=std::min(dim(),o.dim());
        for(int i=0;i<r;i++)
            u[i]+=o.u[i];
        return *this;
    }

    auto& operator-=(const s_vector &o)
    {
        auto r=std::min(dim(),o.dim());
        for(int i=0;i<r;i++)
            u[i]-=o.u[i];
        return *this;
    }

    auto& operator*=(R k)
    {
        for(auto &s:u)
            s*=k;
        return *this;
    }

    auto operator+(const s_vector &o) const
    {
        auto v=*this;
        return v+=o;
    }

    auto operator-(const s_vector &o) const
    {
        auto v=*this;
        return v-=o;
    }

    auto operator-() const
    {
        auto v=*this;
        for(auto &s:v.u)
            s=-s;
        return v;
    }

    auto& operator/=(R k)
    {
        for(auto &s:u)
            s/=k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto v=*this;
        return v/=k;
    }

    auto begin()
    {
        return u.begin();
    }

    auto begin() const
    {
        return u.cbegin();
    }

    auto end()
    {
        return u.end();
    }

    auto end() const
    {
        return u.cend();
    }

    template <size_t k>
    auto& get()& {
        return u[k];
    }

    template <size_t k>
    const auto& get() const& {
        return u[k];
    }

    template <size_t k>
    auto&& get() const&& {
        return u[k];
    }

    template <size_t k>
    auto&& get() && {
        return u[k];
    }
};

template<typename R,int n>
auto operator*(const R&k,const s_vector<R,n>& u)
{
    auto v=u;
    return v*=k;
}

namespace std
{
    template<typename R,int n>
    struct tuple_size<s_vector<R, n>> : std::integral_constant<size_t, n>{};
    template<unsigned long long int k,typename R,int n>
    struct tuple_element<k, s_vector<R, n>>
    {
        using type = R;
    };
}

struct m_shape
{
    int n,m;
};

/*
* This is the union of R^(n*m) for all n and m
* @Requirements
* R is a commutative ring. Formally, it is the set of matrices over R
*/
template<typename R>
class d_matrix
{
    std::vector<std::vector<R>> M;
public:
    using base_field=R;
    using base_ring=R;
    inline static int n=0,m=0;
    d_matrix(R k=0,m_shape shape ={n,m}):M(shape.n,std::vector<R>(shape.m))
    {
        for(int i=0;i<std::min(shape.n,shape.m);i++)
            M[i][i]=k;
    }
    d_matrix(std::vector<std::vector<R>> &&_M):M(std::move(_M)){}
    d_matrix(std::vector<std::vector<R>> _M) :M(std::move(_M)) {}
    auto row_dim() const
    {
        return M.size();
    }

    auto col_dim() const
    {
        return M.empty()?0:M[0].size();
    };

    auto& operator[](int k)
    {
        return M[k];
    }

    R tr() const
    {
        int m=col_dim(),n=row_dim();
        R r=0;
        for(int i=0;i<std::min(n,m);i++)
            r+=M[i][i];
        return r;
    }

    d_matrix T() const
    {
        int m=col_dim(),n=row_dim();
        d_matrix P(0,m_shape{m,n});
        for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                P.M[j][i]=M[i][j];
        return P;
    }

    d_matrix H() const
    {
        int m = col_dim(), n = row_dim();
        d_matrix P(0, m_shape{ m,n });
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
            P.M[j][i] = conj(M[i][j]);
        return P;
    }

    const auto& operator[](int k) const
    {
        return M[k];
    }

    auto &operator+=(const d_matrix &O)
    {
        int r1=std::min(row_dim(),O.row_dim()),r2=std::min(col_dim(),O.col_dim());
        for(int i=0;i<r1;i++) for(int j=0;j<r2;j++)
            M[i][j]+=O.M[i][j];
        return *this;
    }

    auto &operator-=(const d_matrix &O)
    {
        int r1=std::min(row_dim(),O.row_dim()),r2=std::min(col_dim(),O.col_dim());
        for(int i=0;i<r1;i++) for(int j=0;j<r2;j++)
                M[i][j]+=O.M[i][j];
        return *this;
    }

    auto operator+(const d_matrix &O) const
    {
        auto N=*this;
        return N+=O;
    }

    auto operator-(const d_matrix &O) const
    {
        auto N=*this;
        return N-=O;
    }

    auto operator-() const
    {
        auto N=*this;
        for(auto &row:N.M) for(auto &s:row)
            s=-s;
        return N;
    }

    auto operator*(const d_matrix &O) const
    {
        int n=row_dim(),m=col_dim(),p=O.col_dim();
        d_matrix N(0,m_shape{n,p});
        for(int i=0;i<n;i++) for(int k=0;k<m;k++) for(int j=0;j<p;j++)
            N.M[i][j]+=M[i][k]*O.M[k][j];
        return N;
    }

    auto &operator*=(const d_matrix &O)
    {
        auto N=(*this)*O;
        M.swap(N.M);
        return *this;
    }

    auto & operator*=(R k)
    {
        for(auto &row:M) for(auto &u:row)
            u*=k;
        return *this;
    }

    d_vector<R>operator*(const d_vector<R> &u) const
    {
        int n=row_dim(),m=col_dim();
        d_vector<R> v(v_shape{n});
        for(int j=0;j<m;j++) for(int i=0;i<n;i++)
            v[i]+=M[i][j]*u[j];
        return v;
    }

    auto &operator/=(R k)
    {
        for(auto &row:M) for(auto &u:row)
            u/=k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto N=*this;
        return N/=k;
    }

    auto& operator/=(const d_matrix& O)
    {
        return *this *= O.inv();
    }

    auto operator/(const d_matrix &O) const
    {
        return (*this) * O.inv();
    }

    auto begin()
    {
        return M.begin();
    }

    auto begin() const
    {
        return M.cbegin();
    }

    auto end()
    {
        return M.end();
    }

    auto end() const
    {
        return M.cend();
    }

    auto row_echelon_form() const
    {
        int n=row_dim(),m=col_dim();
        auto P=*this;
        int s=0;
        for(int i=0;i<n;i++)
        {
            int p=s;
            while(p<n && P.M[p][i]==0)
                p++;
            if(p==n)
                continue;
            std::swap(P.M[p],P.M[s]);
            R w=P.M[s][i];
            for(int j=s+1;j<n;j++)
            {
                R r=P.M[j][i]/w;
                for (int k = i; k < m; k++)
                    P.M[j][k]-=r*P.M[i][k];
            }
        }
        return P;
    }

    int rank() const
    {
        auto E=row_echelon_form();
        int r=0;
        int n=row_dim(),m=col_dim();
        for(int i=0,j=0;i<n&&j<m;j++)
            if(E.M[i][j]!=0)
            {
                r++;
                i++;
            }
        return r;
    }

    int nullity() const
    {
        return row_dim()-rank();
    }

    R det() const
    {
        int n=row_dim(),m=col_dim();
        auto P=*this;
        bool invert=false;
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(p<n && P.M[p][i]==0)
                p++;
            if(p==n)
                return 0;
            if(p!=i)
            {
                std::swap(P.M[p], P.M[i]);
                invert=!invert;
            }
            R w=P.M[i][i];
            for(int j=i+1;j<n;j++)
            {
                R r=P.M[j][i]/w;
                for (int k = i; k < m; k++)
                    P.M[j][k]-=r*P.M[i][k];
            }
        }
        R d=1;
        for(int i=0;i<n;i++)
            d*=P.M[i][i];
        return invert?-d:d;
    }


    d_matrix inv() const
    {
        int n=row_dim(),m=col_dim();
        d_matrix P=*this,Q(1,m_shape{n,m});
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(p<n && P.M[p][i]==R(0))
                p++;
            if(p==n)
                continue;
            std::swap(P.M[p], P.M[i]);
            std::swap(Q.M[p],Q.M[i]);
            R w=P.M[i][i];
            for(int j=i+1;j<n;j++)
            {
                R r=P.M[j][i]/w;
                for (int k = 0; k < m; k++)
                {
                    P.M[j][k] -= r*P.M[i][k];
                    Q.M[j][k] -= r*Q.M[i][k];
                }
            }
        }
        for(int i=n-1;i>=0;i--)
        {
            R w=P.M[i][i];
            for(int j=0;j<n;j++)
                Q.M[i][j]/=w;
            for(int k=i-1;k>=0;k--)
            {
                R r=P.M[k][i];
                for (int j = 0; j < n; j++)
                    Q.M[k][j] -= r*Q.M[i][j];
            }
        }
        return Q;
    }

    d_vector<R> solve(d_vector<R> A)
    {
        int n=col_dim(),m=row_dim();
        d_matrix P=*this;
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(p<n && P.M[p][i]==R(0))
                p++;
            if(p==n)
                continue;
            std::swap(P.M[p], P.M[i]);
            std::swap(A[p],A[i]);
            R w=P.M[i][i];
            for(int j=i+1;j<n;j++)
            {
                if(w==R(0))
                    continue;
                R r=P.M[j][i]/w;
                for (int k = 0; k < m; k++)
                    P.M[j][k] -= r*P.M[i][k];
                A[j]-=r*A[i];
            }
        }
        for(int i=n-1;i>=0;i--)
        {
            R w=P.M[i][i];
            if(w==R(0))
                continue;
            A[i]/=w;
            for(int k=i-1;k>=0;k--)
            {
                R r=P.M[k][i];
                A[k] -= r*A[i];
            }
        }
        return A;
    }
};

/*
 * It is the union of L(R^a,R^b) over all a and b where L(E,F) is the set of matrices acting on E with values over F
 * - Addition between 2 matrices are defined with respect to the first matrix's shape
 * - Multiplication between 2 matrices is defined if the shapes are compatible
 * - Multiplication between a matrix and a vector is defined if the shapes are compatible
 * - for all (a,b), the set of all vectors of shape (a,b) is an R-vector space
 * - for all a, the set of all matrices of shape (a,a) is an associative algebra over R
 * */


template<typename R>
d_matrix<R> operator*(const R&a,const d_matrix<R> &M)
{
    auto N=M;
    return N*=a;
}

/*
 * It is an element of the vector space L(R^n,R^m)
 * @Requirements
 * R is a commutative ring
 * @Notes
 * - Multiplication between 2 matrices is defined if the shapes are compatible
 * - Multiplication between a matrix and a vector is defined if the shapes are compatible
 * - It is an associative algebra if n=m
 * */

template<typename R,int n,int m>
class s_matrix
{
    std::array<std::array<R,m>,n> M;
public:
    using base_field=R;
    using base_ring=R;
    s_matrix(R k=0)
    {
        for(int i=0;i<n;i++) for(int j=0;j<m;j++)
            M[i][j]=i==j?k:0;
    }
    s_matrix(std::array<std::array<R,m>,n> _M):M(std::move(_M)){}

    inline static constexpr int row_dim()
    {
        return n;
    }

    s_matrix(const std::vector<std::array<R,n>> &_M)
    {
        int counter=0;
        for(int i=0;i<n;i++) for(int j=0;j<m;j++)
            M[i][j]=_M[i][j];
    }

    inline static constexpr int col_dim()
    {
        return m;
    };

    auto& operator[](int k)
    {
        return M[k];
    }

    const auto& operator[](int k) const
    {
        return M[k];
    }

    R tr() const
    {
        R r=0;
        for(int i=0;i<std::min(n,m);i++)
            r+=M[i][i];
        return r;
    }

    s_matrix<R,m,n> T() const
    {
        s_matrix<R,m,n> P;
        for(int i=0;i<n;i++) for(int j=0;j<m;j++)
            P.M[j][i]=M[i][j];
        return P;
    }

    s_matrix<R, m, n> H() const
    {
        s_matrix<R, m, n> P;
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++)
            P.M[j][i] = conj(M[i][j]);
        return P;
    }

    auto &operator+=(const s_matrix &O)
    {
        for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                M[i][j]+=O.M[i][j];
        return *this;
    }

    auto &operator-=(const s_matrix &O)
    {
        for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                M[i][j]+=O.M[i][j];
        return *this;
    }

    auto operator+(const s_matrix &O) const
    {
        auto N=*this;
        return N+=O;
    }

    auto operator-(const s_matrix &O) const
    {
        auto N=*this;
        return N-=O;
    }

    auto operator-() const
    {
        auto N=*this;
        for(auto &row:N.M) for(auto &s:row)
                s=-s;
        return N;
    }

    template<int p>
    s_matrix<R,n,p> operator*(const s_matrix<R,p,m> &O) const
    {
        s_matrix<R,n,p> N(0);
        for(int i=0;i<n;i++) for(int k=0;k<m;k++) for(int j=0;j<p;j++)
                    N.M[i][j]+=M[i][k]*O.M[k][j];
        return N;
    }

    auto &operator*=(const s_matrix &O)
    {
        static_assert(n==m);
        auto N=(*this)*O;
        M.swap(N.M);
        return *this;
    }

    auto & operator*=(R k)
    {
        for(auto &row:M) for(auto &u:row)
            u*=k;
    }
    s_vector<R,n> operator*(const s_vector<R,m> &u) const
    {
        s_vector<R,n> v;
        for(int j=0;j<m;j++) for(int i=0;i<n;i++)
                v[i]+=M[i][j]*u[j];
        return v;
    }

    auto &operator/=(R k)
    {
        for(auto &row:M) for(auto &u:row)
                u/=k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto N=*this;
        return N/=k;
    }

    auto& operator/=(const s_matrix &O)
    {
        return *this*=O.inv();
    }

    auto operator/(const s_matrix &O) const
    {
        return (*this) * O.inv();
    }

    auto begin()
    {
        return M.begin();
    }

    auto begin() const
    {
        return M.cbegin();
    }

    auto end()
    {
        return M.end();
    }

    auto end() const
    {
        return M.cend();
    }

    auto row_echelon_form() const
    {
        auto P=*this;
        int s=0;
        for(int i=0;i<n;i++)
        {
            int p=s;
            while(p<n && P.M[p][i]==0)
                p++;
            if(p==n)
                continue;
            P.M[p].swap(P.M[s]);
            R w=P.M[s][i];
            for(int j=s+1;j<n;j++)
            {
                R r=P.M[j][i]/w;
                for (int k = i; k < m; k++)
                    P.M[j][k]-=r*P.M[i][k];
            }
            s++;
        }
        return P;
    }

    int rank() const
    {
        auto E=row_echelon_form();
        int r=0;
        for(int i=0,j=0;i<n&&j<m;j++)
            if(E.M[i][j]!=0)
            {
                r++;
                i++;
            }
        return r;
    }

    int nullity() const
    {
        return row_dim()-rank();
    }

    R det() const
    {
        static_assert(n==m);
        auto P=*this;
        bool invert=false;
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(p<n &&P.M[p][i]==0)
                p++;
            if(p==n)
                continue;
            if(p!=i)
            {
                std::swap(P.M[p], P.M[i]);
                invert=!invert;
            }
            R w=P.M[i][i];
            for(int j=i+1;j<n;j++)
            {
                R r=P.M[j][i]/w;
                for (int k = i; k < m; k++)
                    P.M[j][k]-=r*P.M[i][k];
            }
        }
        R d=1;
        for(int i=0;i<n;i++)
            d*=P.M[i][i];
        return invert?-d:d;
    }

    s_matrix inv() const
    {
        static_assert(n==m);
        s_matrix P=*this,Q(1);
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(p<n &&P.M[p][i]==R(0))
                p++;
            if(p==n)
                continue;
            std::swap(P.M[p], P.M[i]);
            std::swap(Q.M[p],Q.M[i]);
            R w=P.M[i][i];
            for(int j=i+1;j<n;j++)
            {
                R r=P.M[j][i]/w;
                for (int k = 0; k < m; k++)
                {
                    P.M[j][k] -= r*P.M[i][k];
                    Q.M[j][k] -= r*Q.M[i][k];
                }
            }
        }
        for(int i=n-1;i>=0;i--)
        {
            R w=P.M[i][i];
            for(int j=0;j<n;j++)
                Q.M[i][j]/=w;
            for(int k=i-1;k>=0;k--)
            {
                R r=P.M[k][i];
                for (int j = 0; j < n; j++)
                    Q.M[k][j] -= r*Q.M[i][j];
            }
        }
        return Q;
    }

    s_vector<R,m> solve(s_vector<R,n> A) const
    {
        static_assert(n==m);
        s_matrix P=*this;
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(p<n && P.M[p][i]==R(0))
                p++;
            if(p==n)
                continue;
            std::swap(P.M[p], P.M[i]);
            std::swap(A[p],A[i]);
            R w=P.M[i][i];
            for(int j=i+1;j<n;j++)
            {
                if(w==R(0))
                    continue;
                R r=P.M[j][i]/w;
                for (int k = 0; k < m; k++)
                    P.M[j][k] -= r*P.M[i][k];
                A[j]-=r*A[i];
            }
        }
        for(int i=n-1;i>=0;i--)
        {
            R w=P.M[i][i];
            if(w==R(0))
                continue;
            A[i]/=w;
            for(int k=i-1;k>=0;k--)
            {
                R r=P.M[k][i];
                A[k] -= r*A[i];
            }
        }
        return A;
    }

};

template<typename R,int n,int m>
s_matrix<R,n,m> operator*(const R&a,const s_matrix<R,n,m> &M)
{
    auto N=M;
    return N*=a;
}

template<typename R>
polynomial<R> faddev_lerrier_characteristic_polynomial(const d_matrix<R>&A)
{
    int n=A.row_dim();
    std::vector<R> S(n + 1);
    S[n] = 1;
    d_matrix<R> C(0,m_shape{n,n});
    for (int i = n - 1; i >= 0; i--)
    {
        for (int j = 0; j < n; j++)
            C[j][j] += S[i + 1];
        C = A * C;
        S[i] = -C.tr() / R(n - i);
    }
    return S;
}

template<typename R,int n>
polynomial<R> faddev_lerrier_characteristic_polynomial(const s_matrix<R,n,n>&A)
{
    std::vector<R> S(n + 1);
    S[n] = 1;
    s_matrix<R,n,n> C;
    for (int i = n - 1; i >= 0; i--)
    {
        for (int j = 0; j < n; j++)
            C[j][j] += S[i + 1];
        C = A * C;
        S[i] = -C.tr() / R(n - i);
    }
    return S;
}

template<typename R>
polynomial<R> interpolation_characteristic_polynomial(d_matrix<R> M)
{
    int n=M.row_dim();
    std::vector<R> X(n+1), Y(n+1);
    for (int i = 0; i <= n; i++)
    {
        X[i] = i;
        Y[i] = M.det();
        for (int j = 0; j < n; j++)
            M[j][j] = M[j][j] - 1;
    }
    return newton_interpolation(X, Y);
}

template<typename R,int n>
polynomial<R> interpolation_characteristic_polynomial(s_matrix<R,n,n> M)
{
    std::vector<R> X(n+1), Y(n+1);
    for (int i = 0; i <= n; i++)
    {
        X[i] = i;
        Y[i] = M.det();
        for (int j = 0; j < n; j++)
            M[j][j] = M[j][j] - 1;
    }
    return newton_interpolation(X, Y);
}

template<int n=-1>
using IE = std::conditional_t<n >= 0, s_vector<real, n>, d_vector<real>>;

template<int n = -1>
using IH = std::conditional_t<n >= 0, s_vector<IC, n>, d_vector<IC>>;

template<int n = -1,int m=n>
using IM_IR = std::conditional_t<n >= 0 && m>=0, s_matrix<real, n,m>, d_matrix<real>>;

template<int n = -1, int m = n>
using IM_IC = std::conditional_t<n >= 0 && m >= 0, s_matrix<IC, n, m>, d_matrix<IC>>;
#endif
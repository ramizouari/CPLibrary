//
// Created by ramizouari on 30/11/2021.
//
#include <vector>
#include <array>

using integer = std::int64_t;
using real = long double;
using natural = std::uint64_t;

template<typename R>
R commutator(R a,R b)
{
    return a*b-b*a;
}

template<typename R>
R pow(R a,unsigned long long n)
{
    if(n==0)
        return 1;
    else if(n==1)
        return a;
    auto s=pow(a,n/2);
    return n%2?s*s*a:s*s;
}

struct v_shape
{
    int n;
};

template<typename R>
class d_vector
{
    std::vector<R> u;
public:
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

    auto cbegin() const
    {
        return u.cbegin();
    }

    auto end()
    {
        return u.end();
    }

    auto cend() const
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

template<typename R,int n>
class s_vector
{
    std::array<R,n> u;
public:
    inline static constexpr int dim()
    {
        return n;
    }

    s_vector():u(n){}

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

    auto cbegin() const
    {
        return u.cbegin();
    }

    auto end()
    {
        return u.end();
    }

    auto cend() const
    {
        return u.cend();
    }
};

template<typename R,int n>
auto operator*(const R&k,const s_vector<R,n>& u)
{
    auto v=u;
    return v*=k;
}

struct m_shape
{
    int n,m;
};

template<typename R>
class d_matrix
{
    std::vector<std::vector<R>> M;
public:
    inline static int n=0,m=0;
    d_matrix(R k=0,m_shape shape ={n,m}):M(shape.n,std::vector<R>(shape.m))
    {
        for(int i=0;i<std::min(n,m);i++)
            M[i][i]=k;
    }
    d_matrix(std::vector<std::vector<R>> _M):M(std::move(_M)){}
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

    auto operator*(const d_matrix &O) const
    {
        int n=row_dim(),m=col_dim(),p=O.col_dim();
        d_matrix N(0,m_shape{n,p});
        for(int i=0;i<n;i++) for(int k=0;k<p;k++) for(int j=0;j<m;j++)
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

    auto begin()
    {
        return M.begin();
    }

    auto cbegin() const
    {
        return M.cbegin();
    }
    auto end()
    {
        return M.end();
    }

    auto cend() const
    {
        return M.cend();
    }
    auto row_echelon_form() const
    {
        int n=row_dim(),m=col_dim();
        auto P=*this;
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(P.M[p][i]==0)
                p++;
            if(p==n)
                continue;
            std::swap(P.M[p],P.M[i]);
            R w=P.M[i][i];
            for(int j=i+1;j<n;j++)
            {
                R r=P.M[j][i]/w;
                for (int k = i; k < m; k++)
                    P.M[j][k]-=r*P.M[i][k];
            }
        }
        return P;
    }

    R det() const
    {
        int n=row_dim(),m=col_dim();
        auto P=*this;
        bool invert=false;
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(P.M[p][i]==0)
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


    d_matrix inv() const
    {
        int n=row_dim(),m=col_dim();
        d_matrix P=*this,Q(1,m_shape{n,m});
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(P.M[p][i]==0)
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
};

template<typename R>
d_matrix<R> operator*(const R&a,const d_matrix<R> &M)
{
    auto N=M;
    return N*=a;
}

template<typename R,int n,int m>
class s_matrix
{
    std::array<std::array<R,m>,n> M;
public:
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

    s_matrix(std::vector<std::array<R,n>> _M)
    {
        int counter=0;
        for(const auto &row:_M)
            M[counter++]=row;
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

    auto begin()
    {
        return M.begin();
    }

    auto cbegin() const
    {
        return M.cbegin();
    }

    auto end()
    {
        return M.end();
    }

    auto cend() const
    {
        return M.cend();
    }

    auto row_echelon_form() const
    {
        auto P=*this;
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(P.M[p][i]==0)
                p++;
            if(p==n)
                continue;
            std::swap(P.M[p],P.M[i]);
            R w=P.M[i][i];
            for(int j=i+1;j<n;j++)
            {
                R r=P.M[j][i]/w;
                for (int k = i; k < m; k++)
                    P.M[j][k]-=r*P.M[i][k];
            }
        }
        return P;
    }

    R det() const
    {
        static_assert(n==m);
        auto P=*this;
        bool invert=false;
        for(int i=0;i<n;i++)
        {
            int p=i;
            while(P.M[p][i]==0)
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
            while(P.M[p][i]==0)
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
};

template<typename R,int n,int m>
s_matrix<R,n,m> operator*(const R&a,const s_matrix<R,n,m> &M)
{
    auto N=M;
    return N*=a;
}

template<typename R>
using matrix=d_matrix<R>;

#include <iostream>
constexpr int n=5,m=n;
int main()
{
    matrix<real>::n=n;
    matrix<real>::m=m;
    matrix<real> A;
    for(auto &row:A) for(auto &u:row)
            std::cin >> u;
    auto B=A.inv();
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
            std::cout << B[i][j] << ' ';
        std::cout << '\n';
    }
    std::cout << A.det();
}
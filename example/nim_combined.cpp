#include <iostream>
//
// Created by ramizouari on 30/11/2021.
//


#include <vector>
#include <array>
//
// Created by ramizouari on 01/12/2021.
//



#include <complex>
#include <functional>
#include <cstdint>

namespace cp
{
    using natural = std::uint64_t;
    using integer = std::int64_t;
    using real = long double;
    using IR=real;
    using IC= std::complex<IR>;
    constexpr real epsilon=1e-6;

    template<typename R>
    R commutator(R a,R b)
    {
        return a*b-b*a;
    }

    template<typename M,typename G=typename M::base_field>
    M conj(M a)
    {
        if constexpr (std::is_same_v<G, IC>)
        {
            if constexpr (std::is_same_v<G, M>)
                return std::conj(a);
            else for (auto& s : a)
                    s = conj<typename std::remove_reference<decltype(s)>::type, G>(s);
        }
        return a;
    }

    template<typename R,typename ...StructureMetaData>
    R pow(R a, long long n,StructureMetaData ... meta_info)
    {
        if(n==0)
            return R(1,meta_info...);
        else if(n==1)
            return a;
        auto s=pow(a,n/2,meta_info...);
        return n%2?s*s*a:s*s;
    }

    template<typename R>
    bool is_zero(const R&a)
    {
        return a==R{};
    }


    inline bool is_zero(const std::complex<long double>&a)
    {
        return std::abs(a) < epsilon;
    }

    inline bool is_zero(const std::complex<double>&a)
    {
        return std::abs(a) < epsilon;
    }

    inline bool is_zero(const real &a)
    {
        return std::abs(a) < epsilon;
    }

    template<typename R>
    R gcd(R a,R b)
    {
        if(a<b)
            std::swap(a,b);
        R q,tmp;
        while(!is_zero(b))
        {
            q=a/b;
            tmp=b;
            b=a-b*q;
            a=tmp;
        }
        return a;
    }

    template<typename R>
    R lcm(const R &a,const R &b)
    {
        return a*b/gcd(a,b);
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
        while(!is_zero(b))
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

    template<typename B>
    B next_gray(B n)
    {
        return n^(n>>1);
    }

    template<typename F,typename R>
    std::pair<integer,integer> floyd_functional_cycle(F && f,R x0)
    {
        /*
         * Find a period
         * */
        R x=x0,y=x;
        integer m=0;
        do
        {
            x=f(x);
            y=f(f(y));
            m++;
        }while(y!=x);
        /*
         * Find offset
         * */
        x=x0,y=x;
        for(int i=0;i<m;i++)
            y=f(y);
        int offset=0;
        while(x!=y)
        {
            x=f(x);
            y=f(y);
            offset++;
        }

        /*
         * Find fundamental period
         * */
        y=f(x);
        integer period=1;
        while(x!=y) {
            y = f(y);
            period++;
        }
        return std::make_pair(period,offset);
    }


    template<typename F,typename R>
    integer functional_period(F &&f, R x)
    {
        /*
        * Find a period
        * */
        R y=x;
        integer m=0;
        do
        {
            x=f(x);
            y=f(f(y));
            m++;
        }while(y!=x);
        return m;
    }
}



//
// Created by ramizouari on 01/12/2021.
//


#include <vector>
#include <map>

namespace cp
{
    /**
* @brief Polynomial class
* @details This is the class of polynomials over a commutative ring R
* @tparam R the type of the coefficients
* @Requirements
* <strong>R</strong> is a commutative ring
* @Notes
* Formally this class is simply R[x]
*/
    template<typename R>
    struct polynomial
    {
        std::vector<R> p;

    public:

        const std::vector<R>& data() const
        {
            return p;
        }
        std::vector<R>& data()
        {
            return p;
        }

        bool operator==(const R& a) const
        {
            if (is_zero(a))
                return p.empty();
            else return degree() == 0 && p.front() == a;
        }

        bool operator==(const polynomial<R> &) const = default;

        polynomial& reduce()
        {
            while(!p.empty() && is_zero(p.back()))
                p.pop_back();
            return *this;
        }
        polynomial(R k):p(1,k)
        {
            reduce();
        }

        polynomial(int k) :p(1, k)
        {
            reduce();
        }

        polynomial()=default;

        polynomial(std::vector<R> _p):p(std::move(_p))
        {
            reduce();
        }

        int degree() const
        {
            return p.size()-1;
        }

        int dim() const
        {
            return p.size();
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
            if(is_zero(a))
                p.clear();
            else for(auto& s:p)
                    s*=a;
            reduce();
            return *this;
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

        /**
        * @details creates a preorder between polynomials based on the degree
        * @Requirements:
        * None
        * @Notes
        * This function is essential for the euclidean algorithm to work
        */
        bool operator<(const polynomial &O) const
        {
            return degree() < O.degree();
        }

        /**
         * @brief Polynomial self division
        * @details Divides the polynomial by a constant and stores the result in itself
        * @Requirements
        * One of the following:
         * <ul>
        * <li> R is an integral ring [2]
        * <li> k is invertible
         * </ul>
        */

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

        /**
         * @brief Euclidean division
        * @details Applies euclidean division between two polynomials
        * @Requirements
        * One of the following
        * <ul>
        * <li> R is a field [1]
        * <li> R is an euclidean domain [2]
        * <li> R is a commutative ring, and the dominant coefficient of <strong>O</strong> is inversible
        * <ul/>
        * @Notes
        * Even that condition [1] is a special case of [2], given that some properties of euclidean division are
        * guaranteed only if <strong>R</strong> is a field, We will seperate the two cases
        */
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
            r.reduce();
            return {q,r};
        }

        //Get the quotient of the euclidean division
        polynomial operator/(const polynomial &O) const
        {
            return euclidean_division(O).first;
        }

        //Get the remainder of the euclidean division
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

        polynomial derivative() const
        {
            if (p.empty())
                return {};
            polynomial D;
            D.p.resize(degree());
            for(int i=1;i<=degree();i++)
                D.p[i-1]=p[i]*R(i);
            return D;
        }

        /**
        * @brief Polynomial Evaluation
        * @details Evaluates the polynomial at a given point
        * Evaluates the polynomial over an associative R-algebra H
        * @Requirements
        * <strong>H</strong> is an associative algebra over <strong>R</strong>
        */
        template<typename H>
        H operator()(H a) const
        {
            H r{};
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

        auto begin() const
        {
            return p.begin();
        }

        auto end() const
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

/**
* @brief The functional identity polynomial
* @details This constant is the generator of all polynomials over R.
* @Notes
* Formally, it is simply the polynomial X:x->x
*/
    template<typename R>
    const polynomial<R> X=polynomial<R>(std::vector<R>{0,1});

/**
* @brief Newton Interpolation
* @details Applies Lagrange Interpolation over points (x,y) using Newton's method
* @Requirements
 * <ol>
* <li> x does not have a duplicate element
* <li> One of the following
 * <ul>
*   <li> R is a field
*   <li> (s-t) is invertible for all elements s,t in x
 *  </ul>
 * </ol>
*/
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
}


//
// Created by ramizouari on 09/12/22.
//




#include <vector>
#include <cstddef>
#include <array>

namespace cp::linalg
{
    struct v_shape
    {
        int n;
    };

/**
 * @brief Dynamic Vector
 * @detail Dynamic Vector is a vector in the mathematical sense,
 * @formal It is the union of R^k for all k, where k is the dimension of the vector.
 * <ul>
 * <li> Addition between 2 vectors are defined with respect to the first vector's shape <br>
 * <li> for all k, the set of all vectors of shape k is a vector space
 * </ul>
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

        bool operator==(const d_vector<R>& other) const
        {
            return u==other.u;
        }
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

/**
 * @brief Static Vector:
 * @tparam R is the base field
 * @tparam n is the dimension of the vector space
 * @details It is a member of an R-vector space E where dim(E)= n
 * @Requirements
 * <strong>R</strong> is a commutative ring. <br>
 * @Formal <strong>E</strong> is an <strong>R</strong>-module, and it is a vector space only if <strong>R</strong> is a field. <br>
 * In fact, the name s_vector is used for consistency with the computer science's name.
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

        bool operator==(const s_vector&) const = default;

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
}

namespace std
{
    template<typename R,int n>
    struct tuple_size<cp::linalg::s_vector<R, n>> : std::integral_constant<size_t, n>{};
    template<size_t k,typename R,int n>
    struct tuple_element<k, cp::linalg::s_vector<R, n>>
    {
        using type = R;
    };
}




namespace cp::linalg
{
    struct m_shape
    {
        int n,m;
    };

/**
 * @brief Matrix:
* @details This is the union of R^(n*m) for all n and m
* @Requirements
* <strong>R</strong> is a commutative ring.
* @formal it is the set of matrices over the commutative ring <strong>R</strong> <br>
 * In fact, It is the union of L(R^a,R^b) over all a and b where L(E,F) is the set of matrices acting on E with values over F
*/
    template<typename R>
    struct d_matrix
    {
        std::vector<std::vector<R>> M;

    public:
        using base_field=R;
        using base_ring=R;
        inline static int n=0,m=0;
        const auto& data() const
        {
            return M;
        }
        auto &data()
        {
            return M;
        }
        d_matrix(R k=0,m_shape shape ={n,m}):M(shape.n,std::vector<R>(shape.m))
        {
            for(int i=0;i<std::min(shape.n,shape.m);i++)
                M[i][i]=k;
        }
        d_matrix(std::vector<std::vector<R>> &&_M):M(std::move(_M)){}
        d_matrix(const std::vector<std::vector<R>> &_M) :M(_M) {}
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
                    P.M[j][i] = conj<R,R>(M[i][j]);
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

        auto &operator-=(const R &O)
        {
            int r=std::min(row_dim(),col_dim());
            for(int i=0;i<r;i++)
                M[i][i]-=O;
            return *this;
        }

        auto &operator-=(const d_matrix &O)
        {
            int r1=std::min(row_dim(),O.row_dim()),r2=std::min(col_dim(),O.col_dim());
            for(int i=0;i<r1;i++) for(int j=0;j<r2;j++)
                    M[i][j]-=O.M[i][j];
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

        auto operator-(const R &O) const
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
                while(p<n && is_zero(P.M[p][i]))
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
                s++;
            }
            return P;
        }

        int rank() const
        {
            auto E=row_echelon_form();
            int r=0;
            int n=row_dim(),m=col_dim();
            for(int i=0,j=0;i<n&&j<m;j++)
                if(!is_zero(E.M[i][j]))
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
                while(p<n && is_zero(P.M[p][i]))
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
                while(p<n && is_zero(P.M[p][i]))
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
                while(p<n && is_zero(P.M[p][i]))
                    p++;
                if(p==n)
                    continue;
                std::swap(P.M[p], P.M[i]);
                std::swap(A[p],A[i]);
                R w=P.M[i][i];
                for(int j=i+1;j<n;j++)
                {
                    if(is_zero(w))
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
                if(is_zero(w))
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

    template<typename R>
    d_matrix<R> operator*(const R&a,const d_matrix<R> &M)
    {
        auto N=M;
        return N*=a;
    }

/**
 * @details Static Matrix
 * It is an element of the vector space L(R^n,R^m)
 * @tparam R the base commutative ring
 * @tparam n the number of rows
 * @tparam m the number of columns
 * @Requirements R is a commutative ring
 * @Requirements n>=0
 * @Requirements m>=0
 * @Formal
 * <ul>
 * <li> Multiplication between 2 matrices is defined if the shapes are compatible
 * <li> Multiplication between a matrix and a vector is defined if the shapes are compatible
 * <li> It is an associative algebra if n=m
 * </ul>
 * */

    template<typename R,int n,int m>
    class s_matrix
    {
        std::array<std::array<R,m>,n> M;
    public:
        using base_field=R;
        using base_ring=R;
        bool operator==(const s_matrix&) const = default;
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

        s_matrix(const std::vector<std::array<R,m>> &_M)
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
                    M[i][j]-=O.M[i][j];
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
        s_matrix<R,n,p> operator*(const s_matrix<R,m,p> &O) const
        {
            s_matrix<R,n,p> N(0);
            for(int i=0;i<n;i++) for(int k=0;k<m;k++) for(int j=0;j<p;j++)
                        N[i][j]+=M[i][k]*O[k][j];
            return N;
        }

        auto &operator*=(const s_matrix &O)
        {
            static_assert(n==m);
            auto N=(*this)*O;
            M.swap(N.M);
            return *this;
        }

        auto & operator*=(R k) {
            for (auto &row: M)
                for (auto &u: row)
                    u *= k;
            return *this;
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
                while(p<n && is_zero(P.M[p][i]))
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
                if(! is_zero(E.M[i][j]))
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
                while(p<n && is_zero(P.M[p][i]))
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
                while(p<n && is_zero(P.M[p][i]))
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
                    if(is_zero(w))
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

    template<int n=-1>
    using IE = std::conditional_t<n >= 0, s_vector<real, n>, d_vector<real>>;

    template<int n = -1>
    using IH = std::conditional_t<n >= 0, s_vector<IC, n>, d_vector<IC>>;

    template<int n = -1,int m=n>
    using IM_IR = std::conditional_t<n >= 0 && m>=0, s_matrix<real, n,m>, d_matrix<real>>;

    template<int n = -1, int m = n>
    using IM_IC = std::conditional_t<n >= 0 && m >= 0, s_matrix<IC, n, m>, d_matrix<IC>>;
}

//
// Created by ramizouari on 09/12/22.
//




namespace cp::linalg
{
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

    template<typename IK>
    bool annihilable(const d_matrix<IK> &T, const d_vector<IK> &u,int m)
    {
        d_matrix<IK> A(0,m_shape{m+1,(int)u.dim()});
        A[0]=(std::vector<IK>&)u;
        for(int i=1;i<=m;i++)
        {
            auto v=T * A[i - 1];
            std::copy(v.begin(),v.end(),A[i].begin());
        }
        return A.nullity()!=0;
    }

    template<typename IK>
    polynomial<IK> minimal_polynomial(const d_matrix<IK> &T, const d_vector<IK> &u)
    {
        int n=u.dim();
        std::vector<int> D(n+1);
        for(int i=0;i<=n;i++)
            D[i]=i;
        auto d=*std::upper_bound(D.begin(),D.end(),0,[&T,&u](const auto &x,const auto &y){return annihilable(T,u,x) < annihilable(T,u,y);});
        std::vector<int> mapper(d+1);
        std::vector<polynomial<IK>> Z(d+1);
        std::vector<d_vector<IK>> U(d+1);
        U[0]=u;
        for(int i=1;i<=d;i++)
            U[i] = T * U[i - 1];
        for(int i=0;i<=d;i++)
        {
            mapper[i]=i;
            Z[i].p.resize(d+1);
            Z[i][i]=1;
        }
        int r=0,t=0;
        while(r<=d && t<n)
        {
            if(std::all_of(U[mapper[r]].begin(),U[mapper[r]].end(),[](const auto &x){return is_zero(x);}))
                break;
            int s=r;
            while(s<=d && is_zero(U[mapper[s]][t]))
                s++;
            if(s==d+1)
            {
                t++;
                continue;
            }
            std::swap(mapper[s],mapper[r]);
            Z[mapper[r]]/=U[mapper[r]][t];
            for(int j=t+1;j<n;j++)
                U[mapper[r]][j]/=U[mapper[r]][t];
            U[mapper[r]][t]=1;
            for(int i=r+1;i<=d;i++)
            {
                auto w=U[mapper[i]][t];
                for(int j=t;j<n;j++)
                    U[mapper[i]][j]-=w*U[mapper[r]][j];
                U[mapper[i]][t]=0;
                Z[mapper[i]]-=w*Z[mapper[r]];
            }
            r++;
            t++;
        }
        return Z[mapper[d]]/Z[mapper[d]][Z[mapper[d]].degree()];
    }

    template<typename IK>
    polynomial<IK> minimal_polynomial(const d_matrix<IK> &T, const std::vector<d_vector<IK>> &U)
    {
        polynomial<IK> mu=1;
        for(auto &u:U)
            mu=lcm(mu, minimal_polynomial(T,u));
        return mu;
    }

    template<typename IK>
    polynomial<IK> minimal_polynomial(const d_matrix<IK> &T)
    {
        std::vector<d_vector<IK>> E;
        int n=T.row_dim();
        for(int i=0;i<n;i++)
        {
            E.emplace_back(v_shape{n});
            E.back()[i]=1;
        }
        return minimal_polynomial(T,E);
    }
}




//
// Created by ASUS on 30/11/2021.
//


#include <cstdint>
#include <utility>
//
// Created by ramizouari on 28/11/23.
//





namespace cp
{
    template<integer mod>
    struct cyclic
    {
        integer n;
        inline static bool assume_prime=true;
        inline static constexpr integer m = mod;
        constexpr cyclic(integer o=0):n((o+m)%m){}
        bool operator==(int O) const
        {
            return n==(m+O)%m;
        }

        bool operator!=(int O) const
        {
            return n!=(m+O)%m;
        }

        bool operator==(cyclic O) const
        {
            return n==O.n;
        }

        bool operator!=(cyclic O) const
        {
            return n!=O.n;
        }

        cyclic operator-() const
        {
            return cyclic(-n);
        }

        auto& operator+=(const cyclic &O)
        {
            n=(n+O.n)%m;
            return *this;
        }
        auto& operator-=(const cyclic &O)
        {
            n=(n+m-O.n)%m;
            return *this;
        }

        auto& operator*=(const cyclic &O)
        {
            n=(n*O.n)%m;
            return *this;
        }

        auto& operator/=(const cyclic &O)
        {
            return (*this)*=O.inv();
        }

        auto operator*(const cyclic &O) const
        {
            auto w=*this;
            return w*=O;
        }

        auto operator+(const cyclic &O) const
        {
            auto w=*this;
            return w+=O;
        }

        auto operator-(const cyclic &O) const
        {
            auto w=*this;
            return w-=O;
        }

        auto operator/(const cyclic &O) const
        {
            return (*this)*O.inv();
        }

        cyclic pinv() const
        {
            return egcd(n,m).a;
        }

        auto inv() const
        {
            if(assume_prime)
                return pow(*this,m-2);
            else return pinv();
        }

        auto& operator++()
        {
            return *this+=1;
        }

        auto& operator--()
        {
            return *this-=1;
        }

        auto operator++(int)
        {
            cyclic r(n);
            *this += 1;
            return r;
        }

        auto operator--(int)
        {
            cyclic r(n);
            *this -= 1;
            return r;
        }

        explicit operator integer&()
        {
            return n;
        }

        explicit operator const integer&() const
        {
            return n;
        }

        static constexpr integer modulus()
        {
            return m;
        }
    };

    template<integer m>
    auto operator*(integer k,cyclic<m> s)
    {
        return s*k;
    }

    template<integer m>
    auto operator+(integer k,cyclic<m> s)
    {
        return s+k;
    }

    template<integer m>
    auto operator-(integer k,cyclic<m> s)
    {
        return (-s)+k;
    }
}


//
// Created by ramizouari on 28/11/23.
//



namespace cp
{
    inline static constexpr integer dynamic_modulus=-2;
    template<>
    struct cyclic<dynamic_modulus>
    {
        integer m,n;
        bool assume_prime=true;
    public:
        cyclic(integer o=0,integer q=0):m(q),n(m?(o+m)%m:o){}
        bool operator==(integer O) const
        {
            if(!m) return n==O;
            else return n==(m+O)%m;
        }

        bool operator==(cyclic O) const
        {
            return n==O.n;
        }

        cyclic& operator+=(const cyclic &O)
        {
            if(!m) m=O.m;
            n+=O.n;
            if(m)
                n%=m;
            return *this;
        }

        cyclic& operator+=(integer O)
        {
            n=n+O;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator-=(const cyclic &O)
        {
            if(!m)
                m=O.m;
            n+=m-O.n;
            if(m)
                n%=m;
            return *this;
        }

        cyclic& operator-=(integer O)
        {
            n+=m-O;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator*=(const cyclic &O)
        {
            if(!m) m=O.m;
            n*=O.n;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator*=(integer O)
        {
            n*=O;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator=(integer O)
        {
            n=O;
            if(m) n%=m;
            return *this;
        }

        cyclic inv() const
        {
            if(m==1)
                return *this;
            else if(assume_prime)
                return pow(*this,m-2,m);
            else return pinv();
        }

        cyclic& operator/=(const cyclic &O)
        {
            return (*this)*=O.inv();
        }

        cyclic& operator/=(integer O)
        {
            return (*this)*=cyclic(O,m).inv();
        }

        cyclic operator*(const cyclic &O) const
        {
            auto w=*this;
            return w*=O;
        }

        cyclic operator+(const cyclic &O) const
        {
            auto w=*this;
            return w+=O;
        }

        cyclic operator+(integer O) const
        {
            auto w=*this;
            return w+=O;
        }

        cyclic operator-() const
        {
            return {m-n,m};
        }

        cyclic operator-(const cyclic &O) const
        {
            auto w=*this;
            return w-=O;
        }

        cyclic operator-(integer O) const
        {
            auto w=*this;
            return w-=O;
        }

        cyclic operator/(const cyclic &O) const
        {
            return (*this)*O.inv();
        }

        cyclic operator/(integer O) const
        {
            return (*this)*cyclic(O,m).inv();
        }

        cyclic pinv() const
        {
            return {egcd(n,m).a,m};
        }

        cyclic& operator++()
        {
            return *this+=1;
        }

        cyclic& operator--()
        {
            return *this-=1;
        }

        cyclic operator++(int)
        {
            cyclic r(n,m);
            *this += 1;
            return r;
        }

        cyclic operator--(int)
        {
            cyclic r(n,m);
            *this -= 1;
            return r;
        }

        explicit operator integer&()
        {
            return n;
        }

        explicit operator const integer&() const
        {
            return n;
        }

        integer modulus() const
        {
            return m;
        }
    };

}


//
// Created by ramizouari on 28/11/23.
//



namespace cp
{
    inline static constexpr integer static_modulus=-1;
    template<>
    struct cyclic<static_modulus>
    {
        integer n;
    public:
        inline static integer m=1;
        inline static bool assume_prime=true;
        cyclic(integer o=0):n((o+m)%m){}
        bool operator==(integer O) const
        {
            return n==(m+O)%m;
        }

        bool operator!=(integer O) const
        {
            return n!=(m+O)%m;
        }

        bool operator==(cyclic O) const
        {
            return n==O.n;
        }

        bool operator!=(cyclic O) const
        {
            return n!=O.n;
        }

        cyclic& operator+=(const cyclic &O)
        {
            n=(n+O.n)%m;
            return *this;
        }
        cyclic& operator-=(const cyclic &O)
        {
            n=(n+m-O.n)%m;
            return *this;
        }

        cyclic& operator*=(const cyclic &O)
        {
            n=(n*O.n)%m;
            return *this;
        }

        cyclic inv() const
        {
            if(assume_prime)
                return pow(*this,m-2);
            else return pinv();
        }

        cyclic& operator/=(const cyclic &O)
        {
            return (*this)*=O.inv();
        }

        cyclic operator*(const cyclic &O) const
        {
            auto w=*this;
            return w*=O;
        }

        cyclic operator+(const cyclic &O) const
        {
            auto w=*this;
            return w+=O;
        }

        cyclic operator-() const
        {
            return cyclic(m-n);
        }

        cyclic operator-(const cyclic &O) const
        {
            auto w=*this;
            return w-=O;
        }

        cyclic operator/(const cyclic &O) const
        {
            return (*this)*O.inv();
        }

        cyclic pinv() const
        {
            return egcd(n,m).a;
        }

        cyclic& operator++()
        {
            return *this+=1;
        }

        cyclic& operator--()
        {
            return *this-=1;
        }

        cyclic operator++(int)
        {
            cyclic r(n);
            *this += 1;
            return r;
        }

        cyclic operator--(int)
        {
            cyclic r(n);
            *this -= 1;
            return r;
        }

        explicit operator integer&()
        {
            return n;
        }

        explicit operator const integer&() const
        {
            return n;
        }

        static integer modulus()
        {
            return m;
        }
        static void set_modulus(integer _m)
        {
            m=_m;
        }
    };

    using d_cyclic=cyclic<static_modulus>;

}



//
// Created by ASUS on 01/12/2021.
//



//
// Created by ASUS on 01/12/2021.
//



#include <numbers>
//
// Created by ASUS on 01/12/2021.
//


#include <cstdint>
#include <vector>
#include <map>
#include <numeric>
#include <cmath>
#include <stack>
#include <algorithm>
#include <memory>

namespace cp
{
    using couple =std::pair<integer,integer>;

    class abstract_factoriser
    {
    public:
        virtual ~abstract_factoriser()=default;
        virtual std::vector<integer> prime_factors(integer m) const
        {
            std::vector<integer> P;
            auto p=smallest_divisor(m);
            while(m>1)
            {
                if(P.empty()||P.back()!=p)
                    P.push_back(p);
                m/=p;
            }
            return P;
        }
        virtual integer smallest_divisor(integer m) const=0;
        virtual bool is_prime(integer m) const
        {
            return smallest_divisor(m)==m;
        }
        virtual std::vector<couple> prime_decomposition(integer m) const
        {
            std::vector<couple> P;
            auto p=smallest_divisor(m);
            while(m>1)
            {
                int s=0;
                while(m%p==0)
                {
                    m/=p;
                    s++;
                }
                P.emplace_back(p,s);
            }
            return P;
        }
        virtual std::vector<integer> divisors_list(integer m) const
        {
            std::vector<integer> D;
            divisors_list_rec(m, D, prime_factors(m));
            return D;
        }
    protected:
        void divisors_list_rec(integer n,std::vector<integer> &D,const std::vector<integer> &P, int o=0) const
        {
            auto r=P.size();
            for(int i=o;i<r;i++) if(n%P[i]==0)
                    divisors_list_rec(n/P[i],D,P,i);
            D.push_back(n);

        }

    };

    struct default_factoriser_t
    {
        inline static std::shared_ptr<abstract_factoriser> default_factoriser;
    };


    class factoriser : public abstract_factoriser
    {
        int n;
        std::vector<integer> p_list,smallest_d;
        std::vector<std::vector<integer>> p_factors;


    public:
        factoriser(int _n):n(_n),smallest_d(n+1),p_factors(n+1)
        {
            p_list.reserve(n/log(n));
            std::vector<bool> is_prime(n+1,true);
            for(integer i=2;i<=n;i++) if(is_prime[i])
                {
                    p_list.push_back(i);
                    smallest_d[i]=i;
                    p_factors[i]={i};
                    for(integer j=2*i;j<=n;j+=i)
                    {
                        if(is_prime[j])
                        {
                            smallest_d[j] = i;
                            is_prime[j]=false;
                        }
                        p_factors[j].push_back(i);
                    }
                }
        }

        [[nodiscard]] std::vector<integer> prime_factors(integer m) const override
        {
            if(m<=n)
                return p_factors[m];
            std::vector<integer> result={};
            while(m>1)
            {
                auto p= smallest_divisor(m);
                if(result.empty()||result.back()!=p)
                    result.push_back(p);
                m/=p;
            }
            return result;
        }

        std::vector<std::pair<integer,integer>> prime_decomposition(integer m) const override
        {
            integer r=m;
            std::vector<std::pair<integer,integer>> p_dec;
            for(auto p: prime_factors(m))
            {
                int s=0;
                while(r%p==0)
                {
                    r/=p;
                    s++;
                }
                p_dec.emplace_back(p,s);
            }
            return p_dec;
        }

        std::vector<integer> divisors_list_sorted(integer m) const
        {
            auto D=divisors_list(m);
            std::sort(D.begin(),D.end());
            return D;
        }

        [[nodiscard]] integer smallest_divisor(integer m) const override
        {
            if(m<=n)
                return smallest_d[m];
            integer L=std::ceil(std::sqrt(m));
            for(auto p:p_list)
            {
                if(p>L)
                    break;
                else if(m%p==0)
                    return p;
            }
            return m;
        }

        [[nodiscard]] bool is_prime(integer m) const override
        {
            if(m<n)
                return m>1 && smallest_d[m]==m;
            else
            {
                integer L=std::ceil(std::sqrt(m));
                for(auto p:p_list)
                    if(m%p==0)
                        return false;
                    else if(p>L)
                        break;
                return true;
            }
        }

        integer totient_rec(integer n,const std::vector<integer> &P, integer o=0)
        {
            if(n==0)
                return 0;
            integer S=n;
            for(int i=o;i<P.size();i++)
                S-= totient_rec(n/P[i],P,i+1);
            return S;
        }

        integer totient(integer n)
        {
            integer R=1;
            for(auto [p,m]: prime_decomposition(n))
                R*=pow(p,m-1)*(p-1);
            return R;
        }

        integer totient(integer n,integer m)
        {
            if(n==0)
                return 0;
            auto r=m%n;
            auto P= prime_factors(n);
            return (m/n)*totient(n)+totient_rec(r,P);
        }

        integer carmichael_totient(integer n)
        {
            integer R=1;
            for(auto [p,m]: prime_decomposition(n))
            {
                if(p==2 && m>2)
                    R=std::lcm(R,pow(p,m-2));
                else R=std::lcm(R,pow(p,m-1)*(p-1));
            }
            return R;
        }

        integer divisors_count(integer n)
        {
            integer R=1;
            for(auto [_,m]: prime_decomposition(n))
                R*=(m+1);
            return R;
        }

        integer divisor_function(integer n,integer s)
        {
            if(s==0)
                return divisors_count(n);
            integer R=1;
            for(auto [p,m]: prime_decomposition(n))
                R*=(pow(p,(m+1)*s)-1)/(pow(p,s)-1);
            return R;
        }

        integer divisors_sum(integer n)
        {
            return divisor_function(n,1);
        }

        [[nodiscard]] integer count_primes() const
        {
            return p_list.size();
        }

        [[nodiscard]] const auto& prime_list() const
        {
            return p_list;
        }

        void generate_radicals_rec(std::vector<integer> &R,integer a,integer L,int o=0)
        {
            for(int s=o;s<p_list.size() && a*p_list[s] <= L;s++)
            {
                R.push_back(a*p_list[s]);
                generate_radicals_rec(R,a*p_list[s],L,s+1);
            }
        }

        std::vector<integer> genereate_radicals(integer L)
        {
            std::vector<integer> radicals;
            generate_radicals_rec(radicals,1,L);
            return radicals;
        }

    };

    inline integer chinese_remainder(const std::vector<std::pair<integer,integer>> &S)
    {
        std::stack<std::pair<integer,integer>> Q;
        for(auto s:S)
            Q.push(s);
        while(Q.size() > 1)
        {
            auto [a1,p1]=Q.top();
            Q.pop();
            auto [a2,p2]=Q.top();
            Q.pop();
            auto [k1,k2]=bezout(p1,p2);
            k2*=(a1-a2);
            Q.push({(k2*p2+a2)%(p1*p2),p1*p2});
        }
        return Q.top().first;
    }

    inline integer chinese_remainder(const std::vector<integer>& A,const std::vector<integer>& P)
    {
        std::vector<std::pair<integer,integer>> S;
        int n=A.size(),m=P.size();
        S.reserve(n);
        for(int i=0;i<n;i++)
            S.emplace_back(A[i],P[i]);
        return chinese_remainder(S);
    }
}



#include <algorithm>
#include <optional>
//
// Created by ramizouari on 28/11/23.
//



//
// Created by ramizouari on 01/12/2021.
//




//
// Created by ramizouari on 28/11/23.
//



namespace cp
{
    inline integer totient(integer n,abstract_factoriser &F)
    {
        integer r=1;
        for(auto [p,k]:F.prime_decomposition(n))
            r*=pow(p,k-1)*(p-1);
        return r;
    }

    inline integer carmichael_totient(integer n,abstract_factoriser &F)
    {
        integer r=1;
        for(auto [p,k]:F.prime_decomposition(n))
        {
            if(p==2&&k>=3)
                r=std::lcm(r,pow(2,k-2));
            else
                r = std::lcm(r, pow(p, k - 1) * (p - 1));
        }
        return r;
    }

    inline integer divisors_count(integer n,abstract_factoriser &F)
    {
        integer r=1;
        for(auto [p,k]:F.prime_decomposition(n))
            r*=k+1;
        return r;
    }

    inline integer divisors_sum(integer n,abstract_factoriser &F)
    {
        integer r=1;
        for(auto [p,k]:F.prime_decomposition(n))
            r*=(pow(p,k+1)-1)/(p-1);
        return r;
    }

    inline integer prime_multiplicity(integer n,integer p,abstract_factoriser &F)
    {
        for(auto [q,k]:F.prime_decomposition(n))
            if(q==p)
                return k;
        return 0;
    }

    inline integer prime_multiplicity(integer n,integer p)
    {
        integer r=0;
        while(n%p==0)
        {
            n/=p;
            r++;
        }
        return r;
    }


}


#include <random>

template<cp::integer m>
struct std::hash<cp::cyclic<m>>
{
    inline static std::random_device dev=std::random_device();
    inline static std::mt19937 g=std::mt19937(dev());
    inline static constexpr cp::integer M=1e9+7;
    std::uniform_int_distribution<cp::integer> d=std::uniform_int_distribution<cp::integer>(1,M);
    cp::integer a=d(g),b=d(g);
    public:
    size_t operator()(const cp::cyclic<m> &x) const noexcept
    {
        return (a*static_cast<cp::integer>(x)+b)%M;
    }
};

namespace cp
{
    template<typename cyclic_ring>
    integer discrete_log(cyclic_ring a, cyclic_ring r)
    {
        integer s=std::ceil(std::sqrt(cyclic_ring::m));
        cyclic_ring u=pow(a,s),w=1;
        std::unordered_map<cyclic_ring,integer> mapper;
        for(integer i=0;i<=s;i++,w*=a)
            mapper[r*w]=i;
        w=u;
        for(integer i=1;i<=s;i++,w*=u)
            if(mapper.count(w))
                return i*s-mapper[w];
        return -1;
    }

    inline std::vector<integer> inverse_table(int n,int prime)
    {
        std::vector<integer> I(n + 1);
        I[0] = I[1] = 1;
        for (int i = 2; i <= n; i++)
            I[i] = I[prime % i] *
                   (prime - prime / i) % prime;
        return I;
    }

    inline integer primitive_root_of_unity(integer p,abstract_factoriser &F)
    {
        auto phi=carmichael_totient(p,F);
        auto D=F.divisors_list(phi);
        for(integer k=2;k<p-1;k++) if(std::gcd(k,p)==1)
        {
            bool is_primitive=true;
            for (auto d: D)
                if(d< phi && pow(cyclic<dynamic_modulus>(k,p),d,p)==1)
                {
                    is_primitive=false;
                    break;
                }
            if(is_primitive)
                return k;
        }
        return 0;
    }

    template <integer p>
    integer primitive_root_of_unity(factoriser& F)
    {
        static auto phi = F.totient(p);
        static auto D = F.divisors_list(phi);
        for (integer k = 2; k < p - 1; k++)
        {
            bool is_primitive = true;
            for (auto d : D)
                if (d < phi && pow<d_cyclic>(k, d) == 1)
                {
                    is_primitive = false;
                    break;
                }
            if (is_primitive)
                return k;
        }
        return 0;
    }

    template<integer m>
    integer legendre_symbol(cyclic<m> a)
    {
        integer r;
        if constexpr (m==dynamic_modulus)
            r= (integer) pow(a, (a.modulus() - 1) / 2,a.modulus());
        else
            r= (integer)pow(a, (a.modulus() - 1) / 2);
        if (r > a.modulus() / 2)
            r -= a.modulus();
        return r;
    }

}





//
// Created by ramizouari on 28/11/23.
//





//
// Created by ASUS on 01/12/2021.
//



namespace cp
{
    template<typename T>
    struct binary_operation : public std::enable_shared_from_this<binary_operation<T>>
    {
        using type=T;
        template<typename H0,typename ...H>
        T operator()(const H0&a,const H&... b) const
        {
            if constexpr (sizeof...(b) == 0)
                return a;
            else return reduce(a,this->operator()(b...));
        }
        virtual T reduce(const T& a, const T& b) const = 0;
        virtual T neutral_element() const
        {
            return T{};
        }
    };

    template<typename T>
    struct invertible_operation : public std::enable_shared_from_this<invertible_operation<T>>
    {
        virtual T inv(const T& a) const = 0;
    };

    template<typename T>
    struct multiplicative_inverse_t : public invertible_operation<T>
    {
        T inv(const T& a) const override
        {
            return a.inv();
        }
    };

    template<std::floating_point T>
    struct multiplicative_inverse_t<T> : public invertible_operation<T>
    {
        T inv(const T& a) const override
        {
            return 1/a;
        }
    };

    template<typename T>
    struct multiplicative_inverse_t<std::complex<T>> : public invertible_operation<std::complex<T>>
    {
        std::complex<T> inv(const std::complex<T>& a) const override
        {
            return std::complex<T>(1.)/a;
        }
    };

    template<typename T>
    struct additive_inverse_t : public invertible_operation<T>
    {
        T inv(const T& a) const override
        {
            return -a;
        }
    };

    template<typename T>
    struct monoid_plus_t:public binary_operation<T> {
        T reduce(const T &a, const T &b) const override {
            return a + b;
        }
        inline static T neutral{};
    };

    template<typename T>
    struct plus_t:public monoid_plus_t<T>,public additive_inverse_t<T>
    {
    };

    template<typename T>
    struct multiplies_t:public binary_operation<T>
    {
        T reduce(const T&a,const T&b) const override
        {
            return a*b;
        }

        inline static T neutral=T(1);
        T neutral_element() const override
        {
            return neutral;
        }
    };

    template<typename T>
    struct field_multiplies_t:public multiplies_t<T>,public invertible_operation<T>
    {
    };


    template<typename T>
    class binary_operation_ptr
    {
        std::shared_ptr<binary_operation<T>> op;
    public:
        binary_operation_ptr(std::shared_ptr<binary_operation<T>> value): op(value){}
        template<typename ...H>
        auto operator()(const H&... h) const
        {
            return op->operator()(h...);
        }

        auto reduce(const T& a,const T& b) const
        {
            return op->reduce(a,b);
        }

        auto neutral_element() const
        {
            return op->neutral_element();
        }
    };

    template<typename T>
    class invertible_binary_operation_ptr : public binary_operation_ptr<T>
    {
        std::shared_ptr<invertible_operation<T>> inverter;
    public:
        invertible_binary_operation_ptr(std::shared_ptr<binary_operation<T>> b,
                                        std::shared_ptr<invertible_operation<T>> I): binary_operation_ptr<T>(b),inverter(I){}
        using binary_operation_ptr<T>::operator();
        using binary_operation_ptr<T>::neutral_element;
        auto inv(const T& a) const
        {
            return inverter->inv(a);
        }
    };
}


//
// Created by ASUS on 01/12/2021.
//


#include <vector>
//
// Created by ramizouari on 26/10/23.
//



#include <vector>


//
// Created by ramizouari on 26/10/23.
//



#include <vector>
//
// Created by ramizouari on 26/10/23.
//



namespace cp
{
    inline unsigned int bit_log(unsigned int n)
    {
        unsigned char a=0,b=30,r=0;
        while(a<=b)
        {
            auto c=(a+b)/2;
            if(n>>c)
                a=c+1;
            else
            {
                b=c-1;
                r=c-1;
            }
        }
        if(r && (1<<(r-1))==n)
            return r-1;
        return r;
    }

    inline unsigned int bit_floor(unsigned int n)
    {
        return 1<<bit_log(n);
    }

    inline unsigned int bit_ceil(unsigned int n)
    {
        unsigned r=1;
        while(r<n)
            r<<=1;
        return r;
    }
}



//
// Created by ramizouari on 26/10/23.
//



#include <vector>
namespace cp::data_structures::fixed
{
    template<typename O>
    struct segment_tree
    {
        using R=typename O::type;
        using type=R;
        std::vector<std::vector<R>> S;
        std::vector<R> A;
        int n,h;
        segment_tree(const std::vector<R> &_A):A(_A)
        {
            n=bit_ceil(A.size());
            A.resize(n,O::neutral);
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
                S[m][i]=F(S[m+1][2*i],S[m+1][2*i+1]);
                m--;
                i/=2;
            }
        }

        R query(int l,int r)
        {
            return query(std::max(l,0),std::min(r,n),0,n,0);
        }
    private:
        inline static O F=O();
        void build()
        {
            for(int i=0;i<n;i++)
                S.back()[i]=A[i];
            for(int i=h-2;i>=0;i--) for(int k=0;k<(1<<i);k++)
                    S[i][k]=F(S[i+1][2*k],S[i+1][2*k+1]);
        }
        R query(int l,int r,int a,int b,int depth)
        {
            if(l>=r)
                return O::neutral;
            if(l==a && r==b)
                return S[depth][l>>(h-1-depth)];
            int mid=(a+b)/2;
            if(mid>r)
                return query(l,r,a,mid,depth+1);
            else if(mid<l)
                return query(l,r,mid,b,depth+1);
            else
                return F(query(l,mid,a,mid,depth+1),query(mid,r,mid,b,depth+1));
        }
    };
}


//
// Created by ramizouari on 26/10/23.
//



#include <vector>
namespace cp::data_structures::fixed
{
    template<typename T,typename O>
    struct fenwick_tree {
        int n;
        std::vector<T> bit;
        inline static O F = O();

        fenwick_tree(int _n):n(_n),bit(n,O::neutral){}
        fenwick_tree(const std::vector<T> &X) : fenwick_tree(X.size())
        {
            for(int i=0;i<n;i++)
                update(i,X[i]);
        }
        T sum(int x) {
            if(x<0)
                return O::neutral;
            T ret = O::neutral;
            for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
                ret = F(ret,bit[i]);
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
                bit[i] = F(bit[i], delta);
        }

        void update(int x, T delta) {
            add(x,F(F.inv(sum(x,x)),delta));
        }
    };
}





//
// Created by ramizouari on 27/11/23.
//



//
// Created by ramizouari on 27/11/23.
//



#include <array>
#include <utility>
#include <vector>
#include <numeric>
namespace cp::linalg
{
    inline constexpr std::size_t dynamic_extent = -1;
    template<typename R,std::size_t Rank>
    struct tensor_subview;
    template<typename R,std::size_t Rank>
    struct tensor_view
    {
        virtual R& at(std::array<std::size_t,Rank> indexes) = 0;
        virtual const R& at(std::array<std::size_t,Rank> indexes) const = 0;
        virtual ~tensor_view()= default;
        template<typename ...Args>
        R& at(Args... args)
        {
            return at(std::array<std::size_t,Rank>{args...});
        }
        template<typename ...Args>
        const R& at(Args... args) const
        {
            return at(std::array<std::size_t,Rank>{args...});
        }

        template<typename ...Args>
        R& operator()(Args... args)
        {
            return at(std::array<std::size_t,Rank>{static_cast<std::size_t>(args)...});
        }

        template<typename ...Args>
        const R& operator()(Args... args) const
        {
            return at(std::array<std::size_t,Rank>{static_cast<std::size_t>(args)...});
        }


        const R& operator()(std::array<std::size_t,Rank> args) const
        {
            return at(std::move(args));
        }

        R& operator()(std::array<std::size_t,Rank> args)
        {
            return at(std::move(args));
        }

        virtual std::array<std::size_t,Rank> shape() const = 0;
        static constexpr std::size_t rank()
        {
            return Rank;
        }
        virtual std::size_t size() const
        {
            auto s=shape();
            return std::accumulate(s.begin(),s.end(),1,std::multiplies<>());
        }

        struct iterator
        {
            tensor_view<R,Rank> &src;
            std::array<std::size_t,Rank> indexes;
            bool is_end=false;
            iterator(tensor_view<R,Rank> &src,std::array<std::size_t,Rank> indexes,bool is_end=false):src(src),indexes(indexes),is_end(is_end){}
            iterator& operator++()
            {
                auto shape=src.shape();
                for(int i=Rank-1;i>=0;i--)
                {
                    if(indexes[i]+1<shape[i])
                    {
                        indexes[i]++;
                        return *this;
                    }
                    else
                        indexes[i]=0;
                }
                is_end=true;
                return *this;
            }
            iterator operator++(int)
            {
                iterator tmp=*this;
                ++(*this);
                return tmp;
            }
            bool operator==(const iterator& rhs) const
            {
                return is_end==rhs.is_end && indexes==rhs.indexes;
            }

            bool operator!=(const iterator& rhs) const
            {
                return !(*this==rhs);
            }

            R& operator*()
            {
                return src.at(indexes);
            }

            R* operator->()
            {
                return &src.at(indexes);
            }
        };
        iterator begin()
        {
            return iterator(*this,std::array<std::size_t,Rank>{},false);
        }
        iterator end()
        {
            auto shape=this->shape();
            return iterator(*this,std::array<std::size_t,Rank>{},true);
        }
        virtual tensor_subview<R,Rank> slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end);
        virtual tensor_subview<R,Rank> slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end,std::array<std::size_t,Rank> step);
    };


    template <typename R,std::size_t Rank>
    struct tensor_subview: public tensor_view<R,Rank>
    {
        tensor_view<R,Rank> &src;
        std::array<std::size_t,Rank> m_start,m_end,m_step;
        tensor_subview(tensor_view<R,Rank> &src,std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end):src(src),m_start(start),m_end(end)
        {
            std::fill(m_step.begin(),m_step.end(),1);
        }
        tensor_subview(tensor_view<R,Rank> &src,std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end,
                       std::array<std::size_t,Rank> step):src(src),m_start(start),m_end(end),m_step(step){}

        R& at(std::array<std::size_t,Rank> indexes) override
        {
            std::array<std::size_t,Rank> new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i]*m_step[i] + m_start[i];
            return src.at(new_indexes);
        }
        const R& at(std::array<std::size_t,Rank> indexes) const override
        {
            std::array<std::size_t,Rank> new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i]*m_step[i] + m_start[i];
            return src.at(new_indexes);
        }
        std::array<std::size_t,Rank> shape() const override
        {
            std::array<std::size_t,Rank> new_shape;
            for(int i=0;i<Rank;i++)
                new_shape[i]=(m_end[i]-m_start[i]+m_step[i]-1)/m_step[i];
            return new_shape;
        }
        tensor_subview<R,Rank> slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end) override
        {
            std::array<std::size_t,Rank> new_start,new_end;
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=this->m_start[i]+start[i]*m_step[i];
                new_end[i]=this->m_start[i]+end[i]*m_step[i];
            }
            return tensor_subview<R,Rank>(src,new_start,new_end,m_step);
        }
        tensor_subview<R,Rank> slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end,std::array<std::size_t,Rank> step) override
        {
            std::array<std::size_t,Rank> new_start,new_end,new_step;
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=this->m_start[i]+start[i]*m_step[i];
                new_end[i]=this->m_start[i]+end[i]*m_step[i];
                new_step[i]=this->m_step[i]*step[i];
            }
            return tensor_subview<R,Rank>(src,new_start,new_end,new_step);
        }
    };

    template<typename R,std::size_t Rank>
    tensor_subview<R,Rank> tensor_view<R,Rank>::slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end)
    {
        return tensor_subview<R,Rank>(*this,start,end);
    }

    template<typename R,std::size_t Rank>
    tensor_subview<R,Rank> tensor_view<R,Rank>::slice(std::array<std::size_t,Rank> start,std::array<std::size_t,Rank> end,std::array<std::size_t,Rank> step)
    {
        return tensor_subview<R,Rank>(*this,start,end,step);
    }

    template<typename R>
    struct tensor_view<R,dynamic_extent>
    {
        virtual ~tensor_view()= default;

        template<typename ...Args>
        R& at(Args... args)
        {
            return at(std::vector<std::size_t>{args...});
        }
        template<typename ...Args>
        const R& at(Args... args) const
        {
            return at(std::vector<std::size_t>{args...});
        }

        template<typename ...Args>
        R& operator()(Args... args)
        {
            return at(args...);
        }

        template<typename ...Args>
        const R& operator()(Args... args) const
        {
            return at(args...);
        }
        virtual R& at(std::vector<std::size_t> indexes) = 0;
        virtual const R& at(std::vector<std::size_t> indexes) const = 0;

        [[nodiscard]] virtual std::vector<std::size_t> shape() const = 0;
        virtual std::size_t rank() const
        {
            return shape().size();
        }
        virtual std::size_t size() const
        {
            auto s=shape();
            return std::accumulate(s.begin(),s.end(),1,std::multiplies<>());
        }

        struct iterator
        {
            tensor_view<R,dynamic_extent> &src;
            std::vector<std::size_t> indexes;
            bool is_end=false;
            iterator(tensor_view<R,dynamic_extent> &src,std::vector<std::size_t> indexes,bool is_end=false):src(src),indexes(std::move(indexes)),is_end(is_end){}
            iterator& operator++()
            {
                for(int i=src.rank()-1;i>=0;i--)
                {
                    if(indexes[i]+1<src.shape()[i])
                    {
                        indexes[i]++;
                        return *this;
                    }
                    else
                        indexes[i]=0;
                }
                is_end=true;
                return *this;
            }
            iterator operator++(int)
            {
                iterator tmp=*this;
                ++(*this);
                return tmp;
            }
            bool operator==(const iterator& rhs) const
            {
                return is_end==rhs.is_end && indexes==rhs.indexes;
            }

            bool operator!=(const iterator& rhs) const
            {
                return !(*this==rhs);
            }

            R& operator*()
            {
                return src.at(indexes);
            }
        };
        iterator begin()
        {
            auto zeros=shape();
            for(auto &x:zeros)
                x=0;
            return iterator(*this,zeros,false);
        }
        iterator end()
        {
            auto zeros=shape();
            for(auto &x:zeros)
                x=0;
            return iterator(*this,zeros,true);
        }
        virtual tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end);
        virtual tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end,std::vector<std::size_t> step);
    };

    template<typename R>
    struct vector_view : public tensor_view<R,1>
    {
        R* m_data;
        std::size_t m_size;
        vector_view(R* _data,std::size_t _size):m_data(_data),m_size(_size){}
        vector_view(std::vector<R> &v):m_data(v.data()),m_size(v.size()){}
        template<std::size_t N>
        vector_view(std::array<R,N> &v):m_data(v.data()),m_size(v.size()){}
        std::size_t size() const override
        {
            return m_size;
        }
        virtual R& at(std::size_t i)
        {
            return m_data[i];
        }
        virtual const R& at(std::size_t i) const
        {
            return m_data[i];
        }
        virtual ~vector_view()= default;
        R& at(std::array<std::size_t,1> indexes) override
        {
            return at(indexes[0]);
        }
        const R& at(std::array<std::size_t,1> indexes) const override
        {
            return at(indexes[0]);
        }
        std::array<std::size_t,1> shape() const override
        {
            return {size()};
        }
        tensor_subview<R,1> slice(std::array<std::size_t,1> start,std::array<std::size_t,1> end) override
        {
            return tensor_subview<R,1>(*this,start,end);
        }

        tensor_subview<R,1> slice(std::array<std::size_t,1> start,std::array<std::size_t,1> end, std::array<std::size_t,1> step) override
        {
            return tensor_subview<R,1>(*this,start,end,step);
        }

        vector_view& operator=(const std::vector<R>& O)
        {
            for(int i=0;i<size();i++)
                at(i)=O.at(i);
            return *this;
        }
    };

    template<typename R>
    struct tensor_subview<R,dynamic_extent>:tensor_view<R,dynamic_extent>
    {
        tensor_view<R,dynamic_extent> &src;
        std::vector<std::size_t> m_start,m_end,m_step;
        tensor_subview(tensor_view<R,dynamic_extent> &src,std::vector<std::size_t> start,std::vector<std::size_t> end):src(src),m_start(start),m_end(end),m_step(src.rank())
        {
            std::fill(m_step.begin(),m_step.end(),1);
        }
        tensor_subview(tensor_view<R,dynamic_extent> &src,std::vector<std::size_t> start,std::vector<std::size_t> end,
                       std::vector<std::size_t> step):src(src),m_start(start),m_end(end),m_step(step){}

        R& at(std::vector<std::size_t> indexes) override
        {
            std::vector<std::size_t> new_indexes(src.rank());
            for(int i=0;i<src.rank();i++)
                new_indexes[i]=indexes[i] + m_start[i];
            return src.at(new_indexes);
        }
        const R& at(std::vector<std::size_t> indexes) const override
        {
            std::vector<std::size_t> new_indexes(src.rank());
            for(int i=0;i<src.rank();i++)
                new_indexes[i]=indexes[i] + m_start[i];
            return src.at(new_indexes);
        }
        std::vector<std::size_t> shape() const override
        {
            std::vector<std::size_t> new_shape(src.rank());
            for(int i=0;i<src.rank();i++)
                new_shape[i]=m_end[i]-m_start[i];
            return new_shape;
        }
        tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end) override
        {
            std::vector<std::size_t> new_start(src.rank()),new_end(src.rank());
            for(int i=0;i<src.rank();i++)
            {
                new_start[i]=this->m_start[i]+start[i];
                new_end[i]=this->m_start[i]+end[i];
            }
            return tensor_subview<R,dynamic_extent>(src,new_start,new_end);
        }

        tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end,std::vector<std::size_t> step) override
        {
            auto Rank=src.rank();
            std::vector<std::size_t> new_start(Rank),new_end(Rank),new_step(Rank);
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=this->m_start[i]+start[i]*m_step[i];
                new_end[i]=this->m_start[i]+end[i]*m_step[i];
                new_step[i]=this->m_step[i]*step[i];
            }
            return tensor_subview<R,dynamic_extent>(src,new_start,new_end,new_step);
        }

    };

    template<typename R,std::size_t Rank>
    struct to_dynamic_view_t : public tensor_view<R,dynamic_extent>
    {
        tensor_view<R,Rank> &src;
        to_dynamic_view_t(tensor_view<R,Rank> &src):src(src){}
        R& at(std::vector<std::size_t> indexes) override
        {
            std::array<std::size_t,Rank> new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        const R& at(std::vector<std::size_t> indexes) const override
        {
            std::array<std::size_t,Rank> new_indexes;
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        std::vector<std::size_t> shape() const override
        {
            return std::vector<std::size_t>(src.shape().begin(),src.shape().end());
        }
        tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end) override
        {
            std::vector<std::size_t> new_start(Rank),new_end(Rank);
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=start[i];
                new_end[i]=end[i];
            }
            return tensor_subview<R,dynamic_extent>(*this,new_start,new_end);
        }

        tensor_subview<R,dynamic_extent> slice(std::vector<std::size_t> start,std::vector<std::size_t> end, std::vector<std::size_t> step) override
        {
            std::vector<std::size_t> new_start(Rank),new_end(Rank),new_step(Rank);
            for(int i=0;i<Rank;i++)
            {
                new_start[i]=start[i];
                new_end[i]=end[i];
                new_step[i]=step[i];
            }
            return tensor_subview<R,dynamic_extent>(*this,new_start,new_end,new_step);
        }
    };

    template<typename R>
    struct to_dynamic_view_t<R,dynamic_extent>: public tensor_view<R,dynamic_extent>{};

    template<typename R,std::size_t Rank>
    struct to_static_view_t : public tensor_view<R,Rank>
    {
        tensor_view<R,dynamic_extent> &src;
        std::array<std::size_t,Rank> _shape;
        to_static_view_t(tensor_view<R,dynamic_extent> &src,std::array<std::size_t,Rank> _shape):src(src),_shape(_shape){}
        to_static_view_t(tensor_view<R,dynamic_extent> &src):src(src)
        {
            auto s=src.shape();
            for(int i=0;i<std::min(Rank,src.rank());i++)
                _shape[i]=s[i];
        }
        R& at(std::array<std::size_t,Rank> indexes) override
        {
            std::vector<std::size_t> new_indexes(Rank);
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        const R& at(std::array<std::size_t,Rank> indexes) const override
        {
            std::vector<std::size_t> new_indexes(Rank);
            for(int i=0;i<Rank;i++)
                new_indexes[i]=indexes[i];
            return src.at(new_indexes);
        }
        std::array<std::size_t,Rank> shape() const override
        {
            return _shape;
        }
    };

    template<typename R>
    tensor_subview<R,dynamic_extent> tensor_view<R,dynamic_extent>::slice(std::vector<std::size_t> start,std::vector<std::size_t> end)
    {
        return tensor_subview<R,dynamic_extent>(*this,start,end);
    }

    template<typename R>
    tensor_subview<R,dynamic_extent> tensor_view<R,dynamic_extent>::slice(std::vector<std::size_t> start,std::vector<std::size_t> end,std::vector<std::size_t> step)
    {
        return tensor_subview<R,dynamic_extent>(*this,start,end,step);
    }

    template<typename R,std::size_t Rank>
    to_dynamic_view_t<R,Rank> to_dynamic_view(tensor_view<R,Rank> &src)
    {
        return to_dynamic_view_t<R,Rank>(src);
    }

    template<typename R,std::size_t Rank>
    to_static_view_t<R,Rank> to_static_view(tensor_view<R,dynamic_extent> &src,std::array<std::size_t,Rank> _shape)
    {
        return to_static_view_t<R,Rank>(src,_shape);
    }

    template<std::size_t Rank,typename R>
    to_static_view_t<R,Rank> to_static_view(tensor_view<R,dynamic_extent> &src)
    {
        return to_static_view_t<R,Rank>(src);
    }
}


#include <complex>
#include <cstdint>
namespace cp::signals
{

    enum class FFTNormalization
    {
        None,
        Sqrt,
        Normalized,
        Auto
    };

    using cp::linalg::dynamic_extent;

    template<typename R>
    struct abstract_fft
    {
        virtual void transform(linalg::tensor_view<R,1> &&v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            transform(v,inverse,normalization);
        }
        virtual void transform(linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const = 0;
        virtual void transform(linalg::tensor_view<R,dynamic_extent> &&v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            transform(v,inverse,normalization);
        }
        virtual void transform(linalg::tensor_view<R,dynamic_extent> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            if(v.rank()>1)
                throw std::invalid_argument("rank must be one");
            auto x=to_static_view<1,R>(v);
            transform(x,inverse,normalization);
        }
        R& operator()(R &x, bool inverse=false, FFTNormalization normalization=FFTNormalization::None) const
        {
            transform(x,inverse,normalization);
            return x;
        }
        virtual ~abstract_fft()= default;
    };

    template<typename R>
    void normalize(linalg::tensor_view<std::complex<R>,1> &v,FFTNormalization normalized)
    {
        R r;
        switch (normalized)
        {
            case FFTNormalization::None:
                r=1;
                break;
            case FFTNormalization::Sqrt:
                r=std::sqrt(v.size());
                break;
            case FFTNormalization::Normalized:
                r=v.size();
                break;
            case FFTNormalization::Auto:
                throw std::invalid_argument("cannot normalize with auto");
        }
        if(normalized!=FFTNormalization::None) for (std::complex<R> & x : v)
            x /= r;
    }

    template<typename R>
    void normalize(linalg::tensor_view<std::complex<R>,1> &&v,FFTNormalization normalized)
    {
        normalize(v,normalized);
    }

    template<typename R>
    void inplace_fft2(cp::linalg::tensor_view<std::complex<R>,1> & a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        if(normalized==FFTNormalization::Auto)
            normalized=FFTNormalization::Sqrt;
        int n = a.size();
        for (int i = 1, j = 0; i < n; i++)
        {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                swap(a(i), a(j));
        }
        for (int len = 2; len <= n; len <<= 1)
        {
            R ang = 2 * std::numbers::pi / len * (inverse ? -1 : 1);
            std::complex<R> wlen = std::polar<R>(1.,ang);
            for (int i = 0; i < n; i += len)
            {
                std::complex<R> w(1);
                for (int j = 0; j < len / 2; j++)
                {
                    std::complex<R> u = a(i+j), v = a(i+j+len/2) * w;
                    a(i+j) = u + v;
                    a(i+j+len/2) = u - v;
                    w *= wlen;
                }
            }
        }
        normalize(a,normalized);
    }

    template<typename R>
    void inplace_fft2(cp::linalg::tensor_view<std::complex<R>,1> && a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        inplace_fft2(a,inverse,normalized);
    }

    template<typename R>
    struct radix2_fft: public abstract_fft<R>
    {
        using abstract_fft<R>::transform;
        void transform(linalg::tensor_view<R,1> &v,bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
        {
            if(v.size()!= bit_ceil(v.size()))
                throw std::invalid_argument("size of vector must be a power of 2");
            inplace_fft2(v,inverse,normalization);
        }
    };

    template<typename R>
    struct mixed_radix_fft;


    template<std::floating_point Real>
    struct mixed_radix_fft<std::complex<Real>> : public cp::signals::abstract_fft<std::complex<Real>>, protected cp::default_factoriser_t
    {
        using R=std::complex<Real>;
        using cp::signals::abstract_fft<R>::transform;
        std::shared_ptr<cp::abstract_factoriser> F;
        mixed_radix_fft(std::shared_ptr<cp::abstract_factoriser> _F=default_factoriser):F(_F){}
        void transform_rec(cp::linalg::tensor_view<R,1> &v, bool inverse=false, cp::signals::FFTNormalization normalization = cp::signals::FFTNormalization::None) const
        {
            auto n=v.size();
            if(n==1)
                return;
            std::uint32_t p = F->smallest_divisor(n);
            auto q=n/p;
            std::vector<cp::linalg::tensor_subview<R,1>> V;
            for(unsigned i=0;i<p;i++)
                V.push_back(v.slice({i},{n},{p}));
            for(auto &v:V)
                transform_rec(v,inverse,normalization);
            R w=std::polar(1.0,2*std::numbers::pi/n);
            R z=std::polar(1.0,2*std::numbers::pi/p);
            if(inverse)
            {
                w = std::conj(w);
                z = std::conj(z);
            }
            R t=1;
            std::vector<R> result(n);
            for(int i=0;i<p;i++,t*=z)
            {
                R h1=1,h2=1;
                for (int j = 0; j < p; j++,h1*=t,h2*=w)
                {
                    R h3=1;
                    for (int k = 0; k < q; k++,h3*=h2)
                        result[i*q+k] += h1 * h3 * V[j](k);
                }
            }
            for(int i=0;i<n;i++)
                v(i)=result[i];
        }

        void transform(cp::linalg::tensor_view<R,1> &v, bool inverse=false, cp::signals::FFTNormalization normalization = cp::signals::FFTNormalization::None) const override
        {
            transform_rec(v,inverse,normalization);
            normalize(v,normalization);
        }
    };

    template<typename R>
    R mod_operator(R x, R y)
    {
        return (x%y+y)%y;
    }

    //Bluestein's algorithm
    template<typename R>
    std::vector<std::complex<R>> general_fft(const cp::linalg::tensor_view<std::complex<R>,1> &a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        using cp::linalg::vector_view;
        unsigned int n = a.size();
        if(bit_ceil(n)==n)
        {
            std::vector<std::complex<R>> b(n);
            for(int i=0;i<n;i++)
                b[i]=a(i);
            inplace_fft2(vector_view(b),inverse);
            return b;
        }
        //m>=2*n
        auto m=bit_ceil(2*n-1);
        std::vector<std::complex<R>> w(m);
        R ang = std::numbers::pi  / n * (inverse ? 1 : -1);
        w[1]=std::polar<R>(1.,ang);
        w[0]=1;
        for(int i=2;i<m;i++)
            w[i]=w[i-1]*w[1];
        std::vector<std::complex<R>> A(m),B(m),W(m);
        for(size_t i=0;i<n;i++)
        {
            auto r=mod_operator<std::int64_t>(i*i,2*n);
            W[i]=w[r];
            A[i]=a(i)*w[2*n-r];
            B[i]=w[r];
        }
        for(size_t i=1;i<n;i++)
            B[m-i]=B[i];
        inplace_fft2(vector_view(A),false,FFTNormalization::None);
        inplace_fft2(vector_view(B),false,FFTNormalization::None);
        for(size_t i=0;i<m;i++)
            A[i]*=B[i];
        inplace_fft2(vector_view(A),true,FFTNormalization::Normalized);
        for(size_t i=0;i<n;i++)
            A[i]*=std::conj(W[i]);
        A.resize(n);
        normalize(vector_view(A),normalized);
        return A;
    }

    template<typename R>
    struct bluestein_fft : public abstract_fft<R>
    {
        using abstract_fft<R>::transform;
        void transform(linalg::tensor_view<R,1> &v,bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
        {
            auto b=general_fft(v,inverse,normalization);
            for(int i=0;i<v.size();i++)
                v(i)=b[i];
        }
    };

    template<typename R>
    struct chirpz_transform
    {
        radix2_fft<R> fft;
        R z;
        template<typename ...Args>
        chirpz_transform(R _z,Args&&... args):z(_z),fft(std::forward<Args>(args)...)
        {
        }

        std::vector<R> transform(const cp::linalg::tensor_view<R,1> &a)
        {
            using linalg::vector_view;
            auto n=a.size();
            auto m=bit_ceil(2*n-1);
            std::vector<R> A(m),B(m),W(m);
            for(size_t i=0;i<n;i++)
            {
                W[i]=pow(z,i*i);
                A[i]=a(i)*W[i];
            }
            for(size_t i=0;i<n;i++)
                B[i]=R(1)/W[i];
            for(size_t i=1;i<n;i++)
                B[m-i]=B[i];
            fft.transform(vector_view(A),false,FFTNormalization::None);
            fft.transform(vector_view(B),false,FFTNormalization::None);
            for(size_t i=0;i<m;i++)
                A[i]*=B[i];
            fft.transform(vector_view(A),true,FFTNormalization::Normalized);
            for(size_t i=0;i<n;i++)
                A[i]/=W[i];
            A.resize(n);
            return A;
        }

        std::vector<R> transform(const cp::linalg::tensor_view<R,1> &&a)
        {
            return transform(a);
        }

    };

}


//
// Created by ramizouari on 26/10/23.
//



#include <vector>
#include <memory>

namespace cp::data_structures::dynamic
{
    template<typename R>
    struct segment_tree
    {
        std::vector<std::vector<R>> S;
        std::vector<R> A;
        int n,h;
        binary_operation_ptr<R> F;
        segment_tree(const std::vector<R> &_A, std::shared_ptr<binary_operation<R>> _F):A(_A),F(_F)
        {
            n=bit_ceil(A.size());
            A.resize(n,F.neutral_element());
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
                S[m][i]=F(S[m+1][2*i],S[m+1][2*i+1]);
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
                    S[i][k]=F(S[i+1][2*k],S[i+1][2*k+1]);
        }
        R query(int l,int r,int a,int b,int depth)
        {
            if(l>=r)
                return F.neutral_element();
            if(l==a && r==b)
                return S[depth][l>>(h-1-depth)];
            int mid=(a+b)/2;
            if(mid>r)
                return query(l,r,a,mid,depth+1);
            else if(mid<l)
                return query(l,r,mid,b,depth+1);
            else
                return F(query(l,mid,a,mid,depth+1),query(mid,r,mid,b,depth+1));
        }
    };
}



namespace cp
{

    /**
 * @brief Karatsuba multiplication
* @details Applies Karatsuba multiplication between two polynomials
* @Requirements
* None
*/

    template<typename R>
    polynomial<R> karatsuba_multiplication(const polynomial<R> &p,const polynomial<R> &q)
    {
        constexpr int L=64;
        if(std::min(p.degree(),q.degree())<=L)
            return p*q;
        polynomial<R> a1,b1,a2,b2;
        int n=p.degree(),m=q.degree(),r=std::max(n,m)+1;
        std::vector<R> &u1=a1.data(),&u2=a2.data(),
                &v1=b1.data(),&v2=b2.data();
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
        auto &c=h.data();
        c.resize(n+m+1);
        for(int i=0;i<=r1.degree();i++)
            c[i]+=r1[i];
        for(int i=0;i<=r2.degree();i++)
            c[i+r/2]+=r2[i];
        for(int i=0;i<=r3.degree();i++)
            c[i+s]+=r3[i];
        return h;
    }

    template<typename R>
    struct polynomial_operation;
    template<typename R>
    struct polynomial_operation<polynomial<R>> : public binary_operation<polynomial<R>>
    {
        std::shared_ptr<binary_operation<std::vector<R>>> m_op;
        std::shared_ptr<invertible_operation<R>> m_inv;
        polynomial_operation(std::shared_ptr<binary_operation<std::vector<R>>> _op,std::shared_ptr<invertible_operation<R>> _inv):m_op(_op),m_inv(_inv){}
        polynomial_operation(std::shared_ptr<binary_operation<std::vector<R>>> _op):m_op(_op),m_inv(nullptr){}
        polynomial<R> reduce(const polynomial<R>&a,const polynomial<R>&b) const override
        {
            return m_op->reduce(a.data(),b.data());
        }

        R inv(const R&a) const
        {
            return m_inv->inv(a);
        }

        polynomial<R> neutral_element() const override
        {
            return polynomial<R>(m_op->neutral_element());
        }
        std::shared_ptr<binary_operation<std::vector<R>>> underlying_operator() const
        {
            return m_op;
        }
        std::shared_ptr<invertible_operation<R>> scalar_inverter() const
        {
            return m_inv;
        }
    };

    template<typename R>
    struct karatsuba_multiplies_t : public polynomial_operation<R>
    {
        R reduce(const R&a,const R&b) const override
        {
            return karatsuba_multiplication(a,b);
        }
    };

    template<typename R>
    struct fast_multiplies_t;

    template<typename R>
    struct fast_multiplies_t<std::vector<R>> : public binary_operation<std::vector<R>>
    {
        signals::radix2_fft<R> fft;
        template<typename ...Args>
        fast_multiplies_t(Args&&...args):fft(std::forward<Args>(args)...){}
        std::vector<R> reduce(const std::vector<R>&a,const std::vector<R>&b) const override
        {
            if(a.size()==0 || b.size()==0)
                return {};
            std::vector<R> A(a),B(b);
            auto r=std::bit_ceil<unsigned>(A.size()+B.size()-1);
            A.resize(r);
            B.resize(r);
            linalg::vector_view U(A),V(B);
            fft.transform(U,false,signals::FFTNormalization::None);
            fft.transform(V,false,signals::FFTNormalization::None);
            for(int i=0;i<r;i++)
                A[i]*=B[i];
            fft.transform(U,true,signals::FFTNormalization::Normalized);
            A.resize(a.size()+b.size()-1);
            return A;
        }
        inline static std::vector<R> neutral=std::vector<R>{1};
        std::vector<R> neutral_element() const override
        {
            return neutral;
        }
    };

    template<std::floating_point R>
    struct fast_multiplies_t<std::vector<R>> : public binary_operation<std::vector<R>>
    {
        signals::radix2_fft<std::complex<R>> fft;
        std::vector<R> reduce(const std::vector<R>&a,const std::vector<R>&b) const override
        {
            if(a.size()==0 || b.size()==0)
                return {};
            auto r=std::bit_ceil<unsigned>(a.size()+b.size()-1);
            std::vector<std::complex<R>> A(r);
            for(int i=0;i<a.size();i++)
            {
                A[i].real(a[i]);
                A[i].imag(b[i]);
            }
            linalg::vector_view U(A);
            fft.transform(U,false,signals::FFTNormalization::None);
            std::vector<std::complex<R>> C(r);
            for(unsigned i=0;i<r;i++)
            {
                auto j=(r-i)&(r-1);
                C[i]=A[i]*A[i]-std::conj(A[j]*A[j]);
                C[i]*=std::complex<R>(0,-0.25);
            }
            linalg::vector_view V(C);
            fft.transform(V,true,signals::FFTNormalization::Normalized);
            std::vector<R> C_real(r);
            for(int i=0;i<r;i++)
                C_real[i]=C[i].real();
            C_real.resize(a.size()+b.size()-1);
            return C_real;
        }
        inline static std::vector<R> neutral=std::vector<R>{1};
        std::vector<R> neutral_element() const override
        {
            return neutral;
        }
    };

    template<typename R>
    struct fast_multiplies_t<polynomial<R>> : public binary_operation<polynomial<R>>
    {
        fast_multiplies_t<std::vector<R>> fft;
        polynomial<R> reduce(const polynomial<R>&a,const polynomial<R>&b) const override
        {
            return polynomial<R>(fft.reduce(a.data(),b.data()));
        }
        inline static polynomial<R> neutral=polynomial<R>{1};
    };

    template<typename R>
    std::vector<R> fast_multiplication(const std::vector<R> &A,const std::vector<R> &B)
    {
        static fast_multiplies_t<std::vector<R>> multiplies;
        return multiplies.reduce(A,B);
    }

    template<typename R>
    polynomial<R> fast_multiplication(const polynomial<R> &A,const polynomial<R> &B)
    {
        static fast_multiplies_t<polynomial<R>> multiplies;
        return multiplies.reduce(A,B);
    }

    template<typename Real=real,cp::integer decompositions=2, cp::integer m>
    std::vector<cp::cyclic<m>> fast_modular_multiplication_real(const std::vector<cp::cyclic<m>> &a,const std::vector<cp::cyclic<m>> &b)
    {
        if(a.size()==0)
            return {};
        using namespace cp;
        std::array<std::vector<Real>,decompositions> A,B;
        integer block=std::ceil(std::pow<Real>(a.front().modulus(),1./decompositions));
        for(int i=0;i<a.size();i++)
        {
            auto z=static_cast<integer>(a[i]);
            for(int j=0;j<decompositions;j++)
            {
                auto [q,r]=std::div(z,block);
                A[j].push_back(r);
                z=q;
            }
        }
        for(int i=0;i<b.size();i++)
        {
            auto z=static_cast<integer>(b[i]);
            for(int j=0;j<decompositions;j++)
            {
                auto [q,r]=std::div(z,block);
                B[j].push_back(r);
                z=q;
            }
        }
        std::array<std::array<std::vector<Real>,decompositions>,decompositions> C;
        for(int i=0;i<decompositions;i++) for(int j=0;j<decompositions;j++)
            C[i][j]= fast_multiplication(A[i],B[j]);
        std::vector<cyclic<m>> R(a.size()+b.size()-1);
        for(int i=0;i<R.size();i++)
        {
            integer x=0;
            integer t1=1;
            for(int j=0;j<decompositions;j++)
            {
                integer t2=t1;
                for(int k=0;k<decompositions;k++)
                {
                    x+=std::llround(C[j][k][i])%m*t2%m;
                    t2*=block;
                }
                t1*=block;
            }
            R[i]=x;
        }
        return R;
    }


    template<typename R>
    std::vector<R> formal_inv_2(const std::vector<R> &A,int m)
    {
        if(m==1)
            return {R(1)/A.front()};
        auto B=A;
        for(int i=1;i<A.size();i+=2)
            B[i]=-B[i];
        auto C= fast_multiplication(A,B);
        std::vector<R> T;
        T.resize(m/2);
        for(int i=0;i<T.size() && 2*i < C.size();i++)
            T[i]=C[2*i];
        auto S=formal_inv_2(T,m/2);
        std::vector<R> Q;
        Q.resize(m);
        for(int i=0;i<m/2;i++)
            Q[2*i]=S[i];
        return fast_multiplication(B, Q);
    }

    template<typename R>
    std::vector<R> formal_inv_2(const std::vector<R> &A,int m,std::shared_ptr<binary_operation<std::vector<R>>> multiplies,std::shared_ptr<invertible_operation<R>> inv)
    {
        if(m==1)
            return {inv->inv(A.front())};
        auto B=A;
        for(int i=1;i<A.size();i+=2)
            B[i]=-B[i];
        auto C= multiplies->reduce(A,B);
        std::vector<R> T;
        T.resize(m/2);
        for(int i=0;i<T.size() && 2*i < C.size();i++)
            T[i]=C[2*i];
        auto S=formal_inv_2(T,m/2,multiplies,inv);
        std::vector<R> Q;
        Q.resize(m);
        for(int i=0;i<m/2;i++)
            Q[2*i]=S[i];
        return multiplies->reduce(B, Q);
    }

    template<typename R>
    polynomial<R> formal_inv_2(const polynomial<R> &A,int m)
    {
        return formal_inv_2(A.data(),m);
    }

    template<typename R>
    polynomial<R> formal_inv_2(const polynomial<R> &A,int m,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        return formal_inv_2(A.data(),m,multiplies->underlying_operator(),multiplies->scalar_inverter());
    }

    template<typename R>
    std::vector<R> formal_inv(const std::vector<R> &A,int m)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m));
        C.resize(m);
        return C;
    }

    template<typename R>
    std::vector<R> formal_inv(const std::vector<R> &A,int m,std::shared_ptr<binary_operation<std::vector<R>>> multiplies,std::shared_ptr<invertible_operation<R>> inv)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m),multiplies,inv);
        C.resize(m);
        return C;
    }

    template<typename R>
    polynomial<R> formal_inv(const polynomial<R> &A,int m)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m));
        C.data().resize(m);
        return C;
    }

    template<typename R>
    polynomial<R> formal_inv(const polynomial<R> &A,int m,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m),multiplies, multiplies->scalar_inverter());
        C.data().resize(m);
        return C;
    }

    template<typename R>
    std::vector<R> fast_division(std::vector<R> A,std::vector<R> Q)
    {
        if(A.size()<Q.size())
            return {};
        int m=A.size()-Q.size()+1;
        std::reverse(A.begin(),A.end());
        std::reverse(Q.begin(),Q.end());
        auto P= fast_multiplication(A, formal_inv(Q,m));
        P.resize(m);
        std::reverse(P.begin(),P.end());
        return P;
    }

    template<typename R>
    std::vector<R> fast_division(std::vector<R> A,std::vector<R> Q,std::shared_ptr<binary_operation<std::vector<R>>> multiplies,std::shared_ptr<invertible_operation<R>> inv)
    {
        if(A.size()<Q.size())
            return {};
        int m=A.size()-Q.size()+1;
        std::reverse(A.begin(),A.end());
        std::reverse(Q.begin(),Q.end());
        auto P= multiplies->reduce(A, formal_inv(Q,m,multiplies,inv));
        P.resize(m);
        std::reverse(P.begin(),P.end());
        return P;
    }

    template<typename R>
    polynomial<R> fast_division(const polynomial<R> &A,const polynomial<R> &B)
    {
        return fast_division(A.data(),B.data());
    }

    template<typename R>
    polynomial<R> fast_division(const polynomial<R> &A,const polynomial<R> &B, std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        return fast_division(A.data(),B.data(), multiplies->underlying_operator(),multiplies->scalar_inverter());
    }


    template<typename R>
    polynomial<R> fast_mod(const polynomial<R>&A,const polynomial<R>& B)
    {
        auto P= fast_division(A,B);
        auto Z= A - fast_multiplication(B,P);
        Z.data().resize(B.degree());
        return Z;
    }

    template<typename R>
    polynomial<R> fast_mod(const polynomial<R>&A, const polynomial<R>& B, std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        auto P= fast_division(A,B,multiplies);
        auto Z= A - multiplies->reduce(B,P);
        Z.data().resize(B.degree());
        return Z;
    }

    template<typename R>
    std::pair<polynomial<R>,polynomial<R>> fast_euclidean_division(const polynomial<R> &A,const polynomial<R>& B)
    {
        auto P= fast_division(A,B);
        auto Q=A- fast_multiplication(B,P);
        Q.data().resize(B.degree());
        return std::make_pair(P,Q);
    }

    template<typename R>
    polynomial<R> fast_gcd(const polynomial<R> &A,const polynomial<R> &B)
    {
        if(B==R{})
            return A;
        return fast_gcd(B,fast_mod(A,B).reduce());
    }

    template<typename R>
    polynomial<R> fast_polynomial_expansion(const std::vector<R> &X)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::fixed::segment_tree<fast_multiplies_t<polynomial<R>>> S(P);
        return S.S[0][0];
    }

    template<typename R>
    polynomial<R> fast_polynomial_expansion(const std::vector<R> &X,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::dynamic::segment_tree S(P,multiplies);
        return S.S[0][0];
    }


    template<typename R>
    std::vector<R> fast_multi_evaluation(const polynomial<R> &A,const std::vector<R> &X)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::fixed::segment_tree<fast_multiplies_t<polynomial<R>>> S(P);
        std::vector<polynomial<R>> Z(1<<(S.h-1));
        Z[0]=fast_mod(A,S.S[0][0]);
        for(int i=1;i<S.h;i++)
            for(int j=(1<<i)-1;j>=0;j--)
                Z[j]=fast_mod(Z[j>>1],S.S[i][j]);
        std::vector<R> Y;
        Y.reserve(n);
        for(int i=0;i<n;i++)
            Y.push_back(Z[i](R{}));
        return Y;
    }

    template<typename R>
    std::vector<R> fast_multi_evaluation(const polynomial<R> &A,const std::vector<R> &X,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::dynamic::segment_tree S(P,std::dynamic_pointer_cast<binary_operation<polynomial<R>>>(multiplies));
        std::vector<polynomial<R>> Z(1<<(S.h-1));
        Z[0]=fast_mod(A,S.S[0][0],multiplies);
        for(int i=1;i<S.h;i++)
            for(int j=(1<<i)-1;j>=0;j--)
                Z[j]=fast_mod(Z[j>>1],S.S[i][j],multiplies);
        std::vector<R> Y;
        Y.reserve(n);
        for(int i=0;i<n;i++)
            Y.push_back(Z[i](R{}));
        return Y;
    }

}

//
// Created by ramizouari on 28/11/23.
//




#include <memory>
#include <set>

namespace cp::signals
{

    template<integer n>
    struct abstract_ntt : public abstract_fft<cyclic<n>>, protected default_factoriser_t
    {
        using R=cyclic<n>;
        using abstract_fft<R>::transform;
        abstract_ntt(std::shared_ptr<abstract_factoriser> _F=default_factoriser):F(_F){}
        mutable integer version=0;
        mutable std::array<std::unordered_map<integer,R>,2> cache;
        mutable std::optional<integer> phi;
        mutable R w1,w2;
        mutable std::shared_ptr<abstract_factoriser> F;

        void build() const
        {
            version=cyclic<n>::modulus();
            if(!F)
                F=default_factoriser;
            phi= carmichael_totient(cyclic<n>::modulus(),*F);
            w1= primitive_root_of_unity(cyclic<n>::modulus(),*F);
            w2=w1.pinv();
        }
        virtual R root_of_unity(integer size,integer m,bool inverse) const
        {
            if(version!=cyclic<n>::modulus())
                build();
            R w;
            if(cache[inverse].count(size))
                w=cache[inverse][size];
            else
            {
                w=inverse?w2:w1;
                auto [q,r]=std::div(*phi,size);
                if(r!=0)
                    throw std::invalid_argument("size must divide phi(m)");
                w=pow(w,q);
                cache[inverse][size]=w;
            }
            return w;
        }
        virtual ~abstract_ntt()= default;
    };


    template<integer m>
    void cyclic_normalize(linalg::tensor_view<cyclic<m>,1> &v,FFTNormalization normalized)
    {
        cyclic<m> r;
        switch (normalized)
        {
            case FFTNormalization::None:
                r=1;
                break;
            case FFTNormalization::Sqrt:
                //r=cp::sqrt(cyclic<m>(v.size()));
                throw std::invalid_argument("not implemented");
                break;
            case FFTNormalization::Normalized:
                r=v.size();
                break;
        }
        r=r.pinv();
        for(auto &x:v)
            x*=r;

    }

    template<integer m>
    void cyclic_normalize(linalg::tensor_view<cyclic<m>,1> &&v,FFTNormalization normalized)
    {
        cyclic_normalize(v,normalized);
    }



    template<integer m>
    void inplace_ntt2(cp::linalg::tensor_view<cyclic<m>,1> & a, cyclic<m> theta, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        int n = a.size();
        for (int i = 1, j = 0; i < n; i++)
        {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                std::swap(a(i), a(j));
        }
        std::vector<cyclic<m>> W;
        W.push_back(theta);
        for (int len = 2; len <= n; len <<= 1)
            W.push_back(W.back()*W.back());
        W.pop_back();
        for (int len = 2,r=W.size()-1; len <= n; len <<= 1,r--)
        {
            auto wlen = W[r];
            for (int i = 0; i < n; i += len)
            {
                cyclic<m> w(1);
                for (int j = 0; j < len / 2; j++)
                {
                    cyclic<m> u = a(i+j), v = a(i+j+len/2) * w;
                    a(i+j) = u + v;
                    a(i+j+len/2) = u - v;
                    w *= wlen;
                }
            }
        }
        cyclic_normalize(a,normalized);
    }

    template<typename Cyclic>
    void inplace_ntt2(cp::linalg::tensor_view<Cyclic,1> && a, Cyclic theta, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        inplace_ntt2(a,theta,normalized);
    }

    template<integer n>
    struct radix2_fft<cyclic<n>> : public abstract_ntt<n>
    {
        using R=cyclic<n>;
        using abstract_fft<cyclic<n>>::transform;
        using abstract_ntt<n>::root_of_unity;
        using abstract_ntt<n>::build;
        using default_factoriser_t::default_factoriser;
        radix2_fft(std::shared_ptr<abstract_factoriser> _F=default_factoriser):abstract_ntt<n>(_F)
        {
        }
        void transform(linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
        {
            auto w= root_of_unity(v.size(),cyclic<n>::modulus(),inverse);
            if(normalization==FFTNormalization::Auto)
                normalization=inverse?FFTNormalization::Normalized:FFTNormalization::None;
            inplace_ntt2(v, w, normalization);
        }
    };

}


//
// Created by ramizouari on 16/10/22.
//




#include <random>

namespace cp
{
    inline bool rabin_miller_primality_test(integer n, integer _a)
    {
        if (n <= 2)
            return n == 2;
        else if (n % 2 == 0)
            return false;
        integer r = n - 1, h = 0;
        while (r % 2 == 0)
        {
            r /= 2;
            h++;
        }
        integer d = 1;
        cyclic<dynamic_modulus> a(_a,n);
        auto u = pow(a, r);
        if (u - 1 == 0)
            return true;
        for (int i = 0; i <= h; i++, u *= u) if (u + 1 == 0)
                return true;
        return false;
    }

    inline bool rabin_miller(integer n, integer iter = 7)
    {
        static std::random_device dev;
        static std::mt19937_64 g(dev());
        if (n == 1)
            return false;
        std::uniform_int_distribution<integer> d(2, n - 1);
        for (int i = 0; i < iter; i++) if (!rabin_miller_primality_test(n, d(g)))
            return false;
        return true;
    }

    inline bool rabin_miller(integer n, const std::vector<integer>& provers)
    {
        if (n == 1)
            return false;
        for (const auto& d : provers) if (!rabin_miller_primality_test(n, d))
            return false;
        return true;
    }

    inline bool fermat_test(integer n, const std::vector<integer>& provers)
    {
        if (n == 1)
            return false;
        for (const auto &p : provers) if (pow<cyclic<dynamic_modulus>>(p, n - 1,n) != 1)
                return false;
        return true;
    }

    inline integer rho_divisor_method(integer n,const polynomial<integer> &_P,integer x0)
    {
        polynomial<cyclic<dynamic_modulus>> P;
        for(int i=0;i<=_P.degree();i++)
            P.p.emplace_back(_P[i],n);
        P.reduce();
        cyclic<dynamic_modulus> x(x0,n),y=x;
        integer d=1;
        do
        {
            x=P(x);
            y=P(P(y));
            d=gcd(static_cast<integer>(y-x),n);
        }while(d==1);
        return d;
    };


    class fast_factoriser : public abstract_factoriser
    {
        int iters;
        polynomial<integer> P;
        integer x0;
    public:
        fast_factoriser(int iters,const polynomial<integer> &P,integer x0):iters(iters),P(P),x0(x0){}
        [[nodiscard]] integer smallest_divisor(integer n) const
        {
            if(rabin_miller(n,iters))
                return n;
            return smallest_divisor(rho_divisor_method(n,P,x0));
        }

        [[nodiscard]] std::vector<std::pair<integer,integer>> prime_decomposition(integer n) const
        {
            std::map<integer,integer> M;
            while(n>1)
            {
                auto d=smallest_divisor(n);
                int s=0;
                while(n%d==0)
                {
                    n/=d;
                    s++;
                }
                M.emplace(d,s);
            }
            return {M.begin(),M.end()};
        }
    };

    class randomized_fast_factoriser : public abstract_factoriser
    {
        int iters;
        std::random_device dev;
        mutable std::mt19937_64 g{dev()};
        int max_degree;
    public:
        randomized_fast_factoriser(int iters,int max_degree):iters(iters),max_degree(max_degree){}
        [[nodiscard]] integer smallest_divisor(integer n) const
        {
            std::uniform_int_distribution<integer> d(0,n-1);
            std::vector<integer> P(max_degree+1);
            std::generate(P.begin(),P.end(),std::bind(d,g));
            integer x0=d(g);
            if(rabin_miller(n,iters))
                return n;
            return smallest_divisor(rho_divisor_method(n,P,x0));
        }

        [[nodiscard]] std::vector<std::pair<integer,integer>> prime_decomposition(integer n)
        {
            std::map<integer,integer> M;
            while(n>1)
            {
                auto d=smallest_divisor(n);
                int s=0;
                while(n%d==0)
                {
                    n/=d;
                    s++;
                }
                M.emplace(d,s);
            }
            return {M.begin(),M.end()};
        }

    };
}



#include <random>
constexpr int m=1;
constexpr cp::integer M=998244353;
using IK=cp::cyclic<M>;

int main(int argc, char**argv)
{
    std::ios::sync_with_stdio(false);
    int n,q;
    std::cin >> n >> q;
    cp::linalg::d_matrix<IK> A(0,cp::linalg::m_shape{n,n});
    std::uniform_int_distribution<cp::integer> d(0,M-1);
    std::mt19937_64 rng(std::random_device{}());

    if(argc==1) for(int i=0;i<n;i++) for(int j=0;j<n;j++)
        std::cin >> static_cast<cp::integer&>(A[i][j]);
    else for(int i=0;i<n;i++) for(int j=0;j<n;j++)
        A[i][j]=d(rng);
    auto F=std::make_shared<cp::randomized_fast_factoriser>(5,3);
    cp::default_factoriser_t::default_factoriser=F;
    std::vector<cp::linalg::d_vector<IK>> V(m);
    for(int i=0;i<m;i++)
    {
        V[i]=cp::linalg::d_vector<IK>(cp::linalg::v_shape{n});
        for(int j=0;j<n;j++)
            V[i][j]=d(rng);
    }
    cp::polynomial<IK> P;
    if(n>50)
        P=cp::linalg::minimal_polynomial(A,V);
    else
        P=cp::linalg::faddev_lerrier_characteristic_polynomial(A);
    std::vector<cp::integer> Q(q);
    for(auto &query:Q)
        std::cin >> query;
    std::vector<std::vector<IK>> batch((q+m-1)/m);
    for(int i=0;i<q;i++)
        batch[i/m].push_back(Q[i]);
    std::vector<IK> R;
    for(auto &b:batch)
    {
        auto batch_results=cp::fast_multi_evaluation(P,b);
        R.insert(R.end(),batch_results.begin(),batch_results.end());
    }
    for(auto &r:R)
    {
        if(n&1)
            r=-r;
        std::cout << static_cast<cp::integer>(r) << '\n';
    }
}

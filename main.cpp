#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <numeric>
#include <cmath>
#include <algorithm>


using integer=std::int64_t;
using couple = std::pair<integer,integer>;

template<typename R,typename ...StructureMetaData>
R pow(R a, long long n,StructureMetaData ... meta_info)
{
    if(n==0)
        return R(1,meta_info...);
    else if(n==1)
        return a;
    auto s=pow(a,n/2);
    return n%2?s*s*a:s*s;
}


template<typename R,typename O>
struct segment_tree
{
    std::vector<std::vector<R>> S,lazy;
    R w1,w2;
    std::vector<R> A;
    int n,h;
    segment_tree(const std::vector<R> &_A,R w1,R w2):A(_A),w1(w1),w2(w2)
    {
        n=std::bit_ceil(A.size());
        A.resize(n,O::neutral);
        int m=n;
        h=0;
        while(m)
        {
            m/=2;
            h++;
        }
        S.resize(h);
        lazy.resize(h);
        for(int i=0;i<h;i++) {
            S[i].resize(1 << i);
            lazy[i].resize(1<<i);
        }
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

    void updateRange(int l,int r)
    {
        updateRange(std::max(l,0),std::min(r,n),0,n,0,1);
    }

    R query(int l,int r)
    {
        return query(std::max(l,0),std::min(r,n),0,n,0);
    }
private:
    inline static constexpr O F=O();
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
        if(l==a && r==b) {
            S[depth][l >> (h - 1 - depth)] = F(S[depth][l >> (h - 1 - depth)],
                                               lazy[depth][l >> (h - 1 - depth)] * formula(w1, n >> depth));
            lazy[depth][l >> (h - 1 - depth)]=O::neutral;
            return S[depth][l >> (h - 1 - depth)];
        }
        int mid=(a+b)/2;
        lazy[depth+1][a>>(h-2-depth)]=F(lazy[depth+1][a>>(h-2-depth)],lazy[depth][a>>(h-1-depth)]);
        lazy[depth+1][mid>>(h-2-depth)]=F(lazy[depth+1][mid>>(h-2-depth)],lazy[depth][a>>(h-1-depth)]);
        if(mid>r)
            return query(l,r,a,mid,depth+1);
        else if(mid<l)
            return query(l,r,mid,b,depth+1);
        else
            return F(query(l,mid,a,mid,depth+1),query(mid,r,mid,b,depth+1));
    }

    void updateRange(int l,int r,int a,int b,int depth,R w)
    {
        if(l>=r)
            return;
        if(l==a && r==b)
        {
            lazy[depth][l>>(h-1-depth)]=F(lazy[depth][l>>(h-1-depth)],w);
            return;
        }
        int mid=(a+b)/2;

        S[depth][l>>(h-1-depth)]=F(S[depth][l>>(h-1-depth)],lazy[depth][l>>(h-1-depth)]*formula(w1,r-l));
        lazy[depth+1][a>>(h-2-depth)]=F(lazy[depth+1][a>>(h-2-depth)],lazy[depth][a>>(h-1-depth)]);
        lazy[depth+1][mid>>(h-2-depth)]=F(lazy[depth+1][mid>>(h-2-depth)],lazy[depth][a>>(h-1-depth)]);
        lazy[depth][a>>(h-1-depth)]=O::neutral;
        if(mid>r)
            return updateRange(l,r,a,mid,depth+1,w);
        else if(mid<l)
            return updateRange(l,r,mid,b,depth+1,w*pow(w1,mid-a));
        else
        {
            updateRange(l,mid,a,mid,depth+1,w);
            updateRange(mid,r,mid,b,depth+1,w*pow(w1,mid-a));
        }
    }

    R formula(R w, int m)
    {
        return (pow(w,m)-1)/(w-1);
    }
};


template<typename T>
struct binary_operation
{
    template<typename H0,typename ...H>
    T operator()(const H0&a,const H&... b) const
    {
        if constexpr (sizeof...(b) == 0)
            return a;
        else return reduce(a,this->operator()(b...));
    }
    virtual T reduce(const T& a, const T& b) const = 0;
};

template<typename T>
struct invertible_operation
{
    virtual T inv(const T& a) const = 0;
};

template<typename T>
struct monoid_plus_t:public binary_operation<T> {
    T reduce(const T &a, const T &b) const override {
        return a + b;
    }
    inline static T neutral{};
};

template<typename T>
struct plus_t:public monoid_plus_t<T>,public invertible_operation<T>
{
    T inv(const T&a) const override
    {
        return -a;
    }
};

template<integer mod>
class cyclic
{
    integer n;
public:
    inline static bool assume_prime=true;
    inline static constexpr integer m = mod;
    cyclic(int o=0):n((o+m)%m){}
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


    auto inv() const
    {
        return pow(*this,m-2);
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
};

template<integer m>
auto operator*(int k,cyclic<m> s)
{
    return s*k;
}

template<integer m>
auto operator+(int k,cyclic<m> s)
{
    return s+k;
}

template<typename cyclic_field>
integer legendre_symbol(cyclic_field a)
{
    integer r= (integer)pow(a, (cyclic_field::m - 1) / 2);
    if (r > cyclic_field::m / 2)
        r -= cyclic_field::m;
    return r;
}

template<typename R>
struct polynomial
{
    std::vector<R> p;

public:
    bool operator==(R a) const
    {
        if (a == R(0))
            return p.empty();
        else return degree() == 0 && p.front() == a;
    }

    bool operator==(const polynomial<R> &) const = default;

    void reduce()
    {
        while(!p.empty() && p.back()==R(0))
            p.pop_back();
    }
    polynomial(R k):p(1,k)
    {
        reduce();
    }

    polynomial(int k = 0) :p(1, k)
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
        if(a==R(0))
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
        H r(0);
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

template<typename R>
struct extension_polynomial_t
{
    polynomial<R> p;
};

/**
 * @brief Dynamic Ring Extension
 * @details It is simply the union of R[x]/q over all polynomials q
 * @Requirements
 * One of the following:
 *  - R is a commutative ring and q is a monic polynomial
 *  - R is a field and q is not zero
 * @Notes
 * If the polynomial is not specified, its value will be by default that of the public static member extension_polynomial
 * @Warning
 * the value of extension_polynomial should be initialized.
 */
template<typename R>
class d_ring_extension
{
    polynomial<R> p;
    polynomial<R> q;
public:
    void reduce()
    {
        if(p.degree()>=q.degree())
            p=p%q;
        p.reduce();
    }

    inline static polynomial<R> extension_polynomial;
    d_ring_extension(R k=0,extension_polynomial_t<R> ext={extension_polynomial}):p(k),q(ext.p)
    {
        reduce();
    }
    d_ring_extension(const std::vector<R> &_p,extension_polynomial_t<R> ext={extension_polynomial}):p(_p),q(ext.p)
    {
        reduce();
    }

    d_ring_extension(extension_polynomial_t<R> ext):q(ext.p){}

    auto& operator+=(const d_ring_extension &O)
    {
        p+=O.p;
        return *this;
    }

    auto& operator-=(const d_ring_extension &O)
    {
        p-=O.p();
        return *this;
    }

    auto operator*(const d_ring_extension &O) const
    {
        auto q=*this;
        return q*=O;
    }

    auto& operator*=(const d_ring_extension &O)
    {
        p*=O.p;
        reduce();
        return *this;
    }

    auto operator+(const d_ring_extension &O) const
    {
        auto q=*this;
        return q+=O;
    }

    auto operator-(const d_ring_extension &O) const
    {
        auto q=*this;
        return q-=O;
    }

    auto inv() const
    {
        auto [a,b,d]= egcd(p,q);
        return a/d[0];
    }

    auto pinv() const
    {
        return egcd(p,q).a;
    }

    auto& operator/=(const d_ring_extension &O)
    {
        return (*this)*=O.inv();
    }

    auto operator/(const d_ring_extension &O) const
    {
        return (*this)*O.inv();
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

    auto& operator[](int k)
    {
        return p[k];
    }

    const auto& operator[](int k) const
    {
        return p[k];
    }
};

template<typename cyclic_field>
cyclic_field sqrt(cyclic_field n)
{
    if (cyclic_field::m == 2)
        return n;
    cyclic_field a = 2;
    while (legendre_symbol(a*a-n) != -1)
        ++a;
    extension_polynomial_t<cyclic_field> q = { polynomial<cyclic_field>({n-a * a,0,1}) };
    d_ring_extension<cyclic_field> phi(std::vector<cyclic_field>{ a,1 }, q);
    return pow(phi, (cyclic_field::m+1)/2,q)[0];
}

constexpr int L=1e5;
constexpr integer M=1e9+9;
using IK=cyclic<M>;

int main()
{
    int T;
    std::cin >>T;
    IK w1,w2;
    auto sqrt5=sqrt<IK>(5);
    w1=(sqrt5+1)/2,w2=(-sqrt5+1)/2;
    std::vector<IK> Z{1,2,3,4};
    segment_tree<IK,plus_t<IK>> S(Z,2,3);
    std::cout << (integer)S.query(0,4) << '\n';
    S.updateRange(0,2);
    std::cout << (integer)S.query(0,4) << '\n';

}
#include <iostream>
#include <chrono>
#include <cstdint>
#include <vector>
#include <map>
#include <random>
#include <variant>
#include <unordered_map>
//
// Created by ASUS on 30/11/2021.
//
#ifndef __MODULAR__ARITHMETIC__
#define __MODULAR__ARITHMETIC__
#include <cstdint>
#include <utility>
//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_FIXED_H
#define CPLIBRARY_FIXED_H
//
// Created by ramizouari on 01/12/2021.
//

#ifndef ACPC_PREPARATION_ABSTRACT_ALGEBRA_H
#define ACPC_PREPARATION_ABSTRACT_ALGEBRA_H
#include <complex>
#include <functional>
#include <cstdint>
#include <concepts>

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


    template<typename R,typename F,typename Id>
    R functional_pow(R a,long long n,const F& f,const Id& identity)
    {
        if(n==0)
            return identity;
        else if(n==1)
            return a;
        auto s=functional_pow(a,n/2,f,identity);
        return n%2?f(f(s,s),a):f(s,s);
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


#endif //ACPC_PREPARATION_ABSTRACT_ALGEBRA_H


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

#endif //CPLIBRARY_FIXED_H
//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_DYNAMIC_H
#define CPLIBRARY_DYNAMIC_H
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

#endif //CPLIBRARY_DYNAMIC_H
//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_STATIC_H
#define CPLIBRARY_STATIC_H
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

#endif //CPLIBRARY_STATIC_H
#endif

constexpr int L=1e6+1;
using integer=std::int64_t;
using extended_integer=std::variant<int,std::monostate>;
constexpr integer M=1e9+7;
using IK=cp::cyclic<M>;
using natural=unsigned long long;
struct RNG
{
    static inline constexpr natural salt=0x9e3779b97f4a7c13;
    std::mt19937_64 rng;
    RNG():rng(std::random_device{}()^salt){}
    RNG(natural seed):rng(seed^salt){}
    natural operator()()
    {
        return rng();
    }
    using result_type=natural;
    static constexpr result_type min()
    {
        return std::mt19937_64::min();
    }
    static constexpr result_type max()
    {
        return std::mt19937_64::max();
    }
};

struct int_hash
{
    RNG rng;
    std::uniform_int_distribution<integer> dist;
    integer a,b;
    inline static constexpr integer M=998244353;
    int_hash():dist(1,M-1),a(dist(rng)),b(dist(rng)){}

    std::size_t operator()(integer x) const
    {
        return (a*x+b)%M;
    }
};

struct pair_int_hash
{
    RNG rng;
    std::uniform_int_distribution<integer> dist;
    integer a,b,c;
    inline static constexpr integer M=998244353;
    pair_int_hash():dist(1,M-1),a(dist(rng)),b(dist(rng)){}

    std::size_t operator()(const std::pair<integer,integer> & x) const
    {
        return (a*x.first+b*x.second+c)%M;
    }
};


extended_integer operator+(const extended_integer &a,const extended_integer &b)
{
    if(a.index()==0 && b.index()==0)
        return std::get<0>(a)+std::get<0>(b);
    else return std::monostate{};
}

integer minimum_moves(integer x, std::unordered_map<integer,integer,int_hash> &cache)
{
    if(x==0)
        return 0;
    if(cache.count(x))
        return cache[x];
    auto s=std::to_string(x);
    integer r=x;
    for(int i=s.size()-1;i>=0;i--) if(s[i]!='0')
        r=std::min(r,minimum_moves(x-(s[i]-'0'),cache)+1);
    return cache[x]=r;
}

integer minimum_moves(integer x)
{
    std::unordered_map<integer,integer,int_hash> cache;
    return minimum_moves(x,cache);
}

integer naive_count(integer n)
{
    int counter=0;
    for(int i=0;i<=n;i++)
    {
        auto s=std::to_string(i);
        bool admissible=true;
        for(int j=1;j<s.size() && admissible;j++)
            if(s[j]==s[j-1])
                admissible=false;
        counter+=admissible;
    }
    return counter;
}

constexpr integer base=10,limit=18;


integer number_of_admissible_numbers(std::string S,int left=0)
{
    integer R=0;
    if(S.empty())
        return 1;
    integer r=S[0]-'0';
    if(r>left)
        r--;
    if(r) for(int j=1;j<S.size();j++) r*=9;
    R+=r;
    if(S[0]-'0'!=left)
        R+=number_of_admissible_numbers(S.substr(1),S[0]-'0');
    return R;
}

integer number_of_admissible_numbers(integer n)
{
    if(n<0)
        return 0;
    auto S=std::to_string(n);
    auto R=number_of_admissible_numbers(S);
    for(int i=0;i<S.size();i++)
        R+=number_of_admissible_numbers(std::string(i,'9'),0);
    return R;
}

int main()
{
    integer a,b;
    std::cin >> a >> b;
    std::cout << number_of_admissible_numbers(b)-number_of_admissible_numbers(a-1) << '\n';
}

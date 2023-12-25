#include <iostream>
#include <chrono>
#include <cstdint>
#include <vector>
#include <map>
#include <random>
#include <variant>
#include <unordered_map>
#include <set>
#include "nt/modular_arithmetic.h"
#include "algebra/permutation.h"
#include "nt/number_theory.h"
#include "nt/functions.h"

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

struct rotation : public cp::abstract_permutation
{
    int n,k;
    rotation(int _n,int k):n(_n),k(k){}
    int transform(int i) const override
    {
        return (i+k)%n;
    }
    int size() const override
    {
        return n;
    }
};

struct double_rotation : public cp::abstract_permutation
{
    int n,m,p,q;
    double_rotation(int _n,int _m,int _p,int _q):n(_n),m(_m),p(_p),q(_q){}
    int transform(int i) const override
    {
        auto [x,y]=std::div(i,m);
        return (x+p)%n*m+(y+q)%m;
    }
    int size() const override
    {
        return n*m;
    }
};

struct grid_element
{
    enum direction : bool
    {
        down,right
    };
    direction d;
    integer x,y;

};

using IK=cp::cyclic<M>;
std::unordered_map<std::pair<integer,integer>,IK,pair_int_hash> cache;
IK count_blocks(integer n,integer m)
{
    if(n<m)
        std::swap(n,m);
    if(cache.count({n,m}))
        return cache[{n,m}];
    if(n<=0 && m<=0)
        return n==0 && m==0;
    IK R=0;
    for(int i=1;i<=n;i++)
      R+=count_blocks(n-i,m);
    for(int i=1;i<=m;i++)
        R+=count_blocks(n-i,m-i);
    for(int i=1;i<=m;i++)
        R+=count_blocks(n,m-i);
    return cache[{n,m}]=R;
}

namespace cp
{
    using Permutation = std::vector<int>;

    void operator*=(Permutation& p, Permutation const& q) {
        Permutation copy = p;
        for (int i = 0; i < p.size(); i++)
            p[i] = copy[q[i]];
    }

    int count_cycles(Permutation p) {
        int cnt = 0;
        for (int i = 0; i < p.size(); i++) {
            if (p[i] != -1) {
                cnt++;
                for (int j = i; p[j] != -1;) {
                    int next = p[j];
                    p[j] = -1;
                    j = next;
                }
            }
        }
        return cnt;
    }

    IK solve(int n, int m,  bool flip_invariance=false) {
        Permutation p(n*m), p1(n*m), p2(n*m), p3(n*m);
        for (int i = 0; i < n*m; i++) {
            p[i] = i;
            p1[i] = (i % n + 1) % n + i / n * n;
            p2[i] = (i / n + 1) % m * n + i % n;
            p3[i] = (m - 1 - i / n) * n + (n - 1 - i % n);
        }

        std::set<Permutation> s;
        for (int i1 = 0; i1 < n; i1++) {
            for (int i2 = 0; i2 < m; i2++) {
                if(flip_invariance)
                {
                    s.insert(p);
                    p*=p3;
                    s.insert(p);
                    p*=p3;
                }
                else
                    s.insert(p);
                p *= p2;
            }
            p *= p1;
        }

        IK sum = 0;
        for (Permutation const& p : s) {
            sum += pow<IK>(2, count_cycles(p));
        }
        return sum / s.size();
    }
}
IK count_torus_colorings(int n,int m,cp::integer K)
{
    IK R=0;
    for(int i=0;i<n;i++) for(int j=0;j<m;j++)
        R+=cp::pow<IK>(K,double_rotation(n,m,i,j).number_of_cycles());
    return R/=n*m;
}

template<std::integral J, std::integral ...I>
std::common_type_t<I...> projection(J a,I ...b)
{
    return a;
}

template<std::integral J,std::integral ...I>
std::common_type_t<J,I...> gcd_closure(J x,I ...i)
{
    if constexpr (sizeof...(i)==0)
        return x;
    else return std::gcd(x,gcd_closure(i...));
}

IK polya_torus_count_naive(cp::integer n, cp::integer m,cp::integer K)
{
    IK R;
    for(int i=0;i<n;i++) for(int j=0;j<m;j++)
    {
        auto a=std::gcd(i,n),b=std::gcd(j,m),c=std::gcd(n/a,m/b);
        auto s=a*b*c;
        R+=cp::pow<IK>(K,s);
    }
    return R/(n*m);
}

IK polya_torus_count(cp::integer n,cp::integer m,cp::integer K,cp::abstract_factoriser &F)
{
    IK R;
    auto A=F.divisors_list(n),B=F.divisors_list(m);
    for(auto d1:A) for(auto d2:B)
        R+=cp::totient(n/d1,F)*cp::totient(m/d2,F)*cp::pow<IK>(K,n*m/std::lcm(n/d1,m/d2));
    return R/(n*m);
}

int main()
{
    auto F=std::make_shared<cp::factoriser>(L);
    int n,m;
    std::cin >> n >> m;
    auto t1=std::chrono::high_resolution_clock::now();
    std::cout << (cp::integer)polya_torus_count(n,m,2,*F) << std::endl;
    auto t2=std::chrono::high_resolution_clock::now();
    std::cout << "Time 1: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << std::endl;}
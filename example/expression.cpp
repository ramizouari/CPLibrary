#include <iostream>
#include <chrono>
#include <cstdint>
#include <vector>
#include <map>
#include <random>
#include <variant>
#include <unordered_map>
#include "nt/modular_arithmetic.h"

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

int main()
{
    int n;
    std::cin >> n;
    std::vector<std::vector<IK>> H(n+1,std::vector<IK>(n+1));
    H[0][0]=1;
    for(int i=1;i<=n;i++) for(int j=1;j<=i;j++)
    {
        for(int k=1;k<=i;k++)
            H[i][j]+=H[i-k][k];
        for(int k=1;k<=j;k++)
            H[i][j]+=H[i][j-k];
        for(int k=1;k<=std::min(i,j);k++)
            H[i][j]+=H[i-k][j-k];
    }
    std::cout << (integer)H[n][n] << '\n';
    std::cout << (integer)count_blocks(n,n) << '\n';
}
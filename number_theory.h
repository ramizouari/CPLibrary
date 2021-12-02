//
// Created by ASUS on 01/12/2021.
//
#ifndef __NT_H__
#define __NT_H__
#include <cstdint>
#include <vector>
#include <map>
#include <numeric>
#include <cmath>

using integer = std::int64_t;
using couple =std::pair<integer,integer>;

class factoriser
{
    int n;
    std::vector<integer> p_list,smallest_d;
    std::vector<std::vector<integer>> p_factors;
    std::vector<std::vector<integer>> d_list;
    std::vector<std::vector<couple>> p_dec;

    void divisors_list_rec(int n,std::vector<integer> &D,const std::vector<integer> &P, int o=0)
    {
        D.push_back(n);
        auto r=P.size();
        for(int i=i;i<=r;i++) if(n%D[i]==0)
                divisors_list_rec(n/P[i],D,P,i);
    }
public:
    factoriser(int _n):n(_n),smallest_d(n+1),p_factors(n+1),d_list(n+1),p_dec(n+1)
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

    const auto& prime_factors(integer m) const
    {
        return p_factors[m];
    }

    const auto& prime_decomposition(integer m)
    {
        integer r=m;
        if(p_dec[m].empty())
        {
            for(auto p: prime_factors(m))
            {
                int s=0;
                while(r%p==0)
                {
                    r/=p;
                    s++;
                }
                p_dec[m].emplace_back(p,s);
            }
        }
        return p_dec[m];
    }

    const auto& divisors_list(integer m)
    {
        if(d_list[m].empty())
            divisors_list_rec(m,d_list[m],p_factors[m]);
        return d_list[m];
    }

    integer smallest_divisor(integer n) const
    {
        return smallest_d[n];
    }

    bool is_prime(int n) const
    {
        return n>1 && smallest_d[n]==n;
    }

    integer totient(integer n)
    {
        integer R=1;
        for(auto [p,m]: prime_decomposition(n))
            R*=pow(p,m-1)*(p-1);
        return R;
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
        if(n==0)
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

};
#endif
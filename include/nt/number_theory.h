//
// Created by ASUS on 01/12/2021.
//
#ifndef __NUMBERTHEORY_H__
#define __NUMBERTHEORY_H__
#include <cstdint>
#include <vector>
#include <map>
#include <numeric>
#include <cmath>
#include <stack>
#include <algorithm>
#include "algebra/abstract_algebra.h"
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
            while(m>1)
            {
                auto p=smallest_divisor(m);
                P.push_back(p);
                m/=p;
                while(m%p==0)
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
            while(m>1)
            {
                auto p=smallest_divisor(m);
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

    template<std::integral I>
    inline I chinese_remainder(const std::vector<std::pair<I,I>> &S)
    {
        std::stack<std::pair<I,I>> Q;
        for(auto s:S)
            Q.push(s);
        while(Q.size() > 1)
        {
            auto [a1,p1]=Q.top();
            Q.pop();
            auto [a2,p2]=Q.top();
            Q.pop();
            auto [k1,k2]=bezout(p1,p2);
            auto x=(k1*a2*p1+k2*a1*p2)%(p1*p2);
            if(x<0)
                x+=p1*p2;
            Q.emplace(x,p1*p2);
        }
        return Q.top().first;
    }

    template<std::integral I>
    I chinese_remainder(const std::vector<I>& A,const std::vector<I>& P)
    {
        std::vector<std::pair<I,I>> S;
        int n=A.size(),m=P.size();
        S.reserve(n);
        for(int i=0;i<n;i++)
            S.emplace_back(A[i],P[i]);
        return chinese_remainder(S);
    }
}


#endif
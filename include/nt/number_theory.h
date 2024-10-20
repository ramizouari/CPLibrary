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
        virtual integer smallest_maximal_prime_power_divisor(integer m) const
        {
            integer p1= smallest_divisor(m) ,p2= smallest_divisor(m/=p1);
            integer q=p1;
            while(p1==p2 && p1 > 1)
            {
                p1=p2;
                q*=p1;
                p2= smallest_divisor(m/=p1);
            }
            return q;
        }

        virtual integer smallest_prime_multiplicity(integer m) const
        {
            if(m==1)
                return 0;
            integer mu=1;
            auto p = smallest_divisor(m);
            while(m > 1 && smallest_divisor(m/=p) == p)
                mu++;
            return mu;
        }
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

        std::vector<integer> divisors_list_sorted(integer m) const
        {
            auto D=divisors_list(m);
            std::sort(D.begin(),D.end());
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

    class light_factoriser : public abstract_factoriser
    {
        int n;
        std::vector<integer> p_list,smallest_d;
    public:
        light_factoriser(int _n):n(_n),smallest_d(n+1)
        {
            p_list.reserve(n/log(n));
            std::vector<bool> is_prime(n+1,true);
            for(integer i=2;i<=n;i++) if(is_prime[i])
            {
                p_list.push_back(i);
                smallest_d[i]=i;
                for(integer j=2*i;j<=n;j+=i) if(is_prime[j])
                {
                    smallest_d[j] = i;
                    is_prime[j]=false;
                }
            }
        }

        [[nodiscard]] std::vector<integer> prime_factors(integer m) const override
        {
            std::vector<integer> result;
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

        [[nodiscard]] integer count_primes() const
        {
            return p_list.size();
        }

        [[nodiscard]] const auto& prime_list() const
        {
            return p_list;
        }
    };


    class fast_factoriser : public abstract_factoriser
    {
        int n;
        std::vector<integer> p_list,smallest_d, smallest_q;
        std::vector<std::uint8_t> q_multiplicity;
    public:
        fast_factoriser(int _n): n(_n), smallest_d(n + 1), smallest_q(n + 1), q_multiplicity(n + 1)
        {
            p_list.reserve(n/log(n));
            std::vector<bool> is_prime(n+1,true);
            for(integer i=2;i<=n;i++) if(is_prime[i])
                {
                    for(integer j=1,k=i;k<=n;k*=i,j++) if(is_prime[k])
                    {
                        q_multiplicity[k]=j;
                        smallest_q[k]=k;
                        for(integer l=2*k;l<=n;l+=k) if(is_prime[l])
                            smallest_q[l]=k;
                    }
                    p_list.push_back(i);
                    smallest_d[i]=i;
                    for(integer j=2*i;j<=n;j+=i) if(is_prime[j])
                    {
                        smallest_d[j] = i;
                        is_prime[j]=false;
                    }
                }
        }

        [[nodiscard]] integer smallest_maximal_prime_power_divisor(integer m) const override
        {
            if(m<=n)
                return smallest_q[m];
            integer L=std::ceil(std::sqrt(m));
            for(auto p:p_list)
            {
                integer q=1;
                if(p>L)
                    break;
                else if(m%p==0)
                {
                    while(m%p==0)
                    {
                        q*=p;
                        m/=p;
                    }
                    return q;
                }
            }
            return m;
        }

        [[nodiscard]] integer smallest_prime_multiplicity(integer m) const override
        {
            if(m<=n)
                return q_multiplicity[smallest_q[m]];
            integer L=std::ceil(std::sqrt(m));
            for(auto p:p_list)
            {
                integer mu=0;
                if(p>L)
                    break;
                else if(m%p==0)
                {
                    while(m%p==0)
                    {
                        mu++;
                        m/=p;
                    }
                    return mu;
                }
            }
            return 1;
        }

        [[nodiscard]] std::vector<integer> prime_factors(integer m) const override
        {
            std::vector<integer> result;
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
            std::vector<std::pair<integer,integer>> p_dec;
            while(m>n)
            {
                auto d= smallest_divisor(m);
                p_dec.emplace_back(d,1);
                while(m%d==0)
                {
                    m/=d;
                    p_dec.back().second++;
                }
            }
            while(m>1)
            {
                auto q=smallest_q[m];
                p_dec.emplace_back(smallest_d[q],q_multiplicity[q]);
                m/=q;
            }
            return p_dec;
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

        [[nodiscard]] integer count_primes() const
        {
            return p_list.size();
        }

        [[nodiscard]] const auto& prime_list() const
        {
            return p_list;
        }
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

        [[nodiscard]] integer count_primes() const
        {
            return p_list.size();
        }

        [[nodiscard]] const auto& prime_list() const
        {
            return p_list;
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
        return Q.empty()? I{} : Q.top().first;
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
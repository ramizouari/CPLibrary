//
// Created by ramizouari on 11/12/23.
//

#ifndef CPLIBRARY_POLY_SPECIAL_POLYNOMIALS_H
#define CPLIBRARY_POLY_SPECIAL_POLYNOMIALS_H

#include <queue>
#include "polynomial/polynomial.h"
#include "nt/number_theory.h"
#include "polynomial/fast_polynomial.h"

namespace cp
{

    template<typename R>
    struct cmp_degree_t
    {
        bool operator()(const polynomial<R>&x,polynomial<R> &y) const
        {
            return x.degree() > y.degree();
        }
    };

    template<typename R>
    polynomial<R> fast_cyclotomic_polynomial(int n,abstract_factoriser &F,std::vector<polynomial<R>> &cache)
    {

        if(n==1)
            return polynomial<R>({-1,1});
        std::priority_queue<polynomial<R>,std::vector<polynomial<R>>, cmp_degree_t<R>> Q;
        for(auto d:F.divisors_list(n)) if(d!=n)
        {
            if(cache[d].data().empty())
                cache[d]=fast_cyclotomic_polynomial(d,F,cache);
            Q.push(cache[d]);
        }
        while(Q.size() > 1)
        {
            auto a=Q.top();
            Q.pop();
            auto b=Q.top();
            Q.pop();
            Q.push(fast_multiplication(a,b));
        }
        polynomial<R> Z;
        Z.data().resize(n+1);
        Z.data().back()=1;
        Z.data().front()=-1;
        return fast_division(Z,Q.top());
    }

    template<typename R>
    polynomial<R> fast_cyclotomic_polynomial(int n,abstract_factoriser &F)
    {
        std::vector<polynomial<R>> cache(n+1);
        return fast_cyclotomic_polynomial(n,F,cache);
    }


    polynomial<cp::integer> integer_cyclotomic_polynomial(int n,abstract_factoriser &F)
    {
        constexpr cp::integer M=998244353;
        using IK=cp::cyclic<M>;
        auto R= fast_cyclotomic_polynomial<IK>(n,F);
        polynomial<cp::integer> P;
        P.data().resize(R.data().size());
        for(int i=0;i<R.data().size();i++) {
            P.data()[i] = static_cast<cp::integer>(R.data()[i]);
            if(P.data()[i]>=M/2)
                P.data()[i]-=M;
        }
        return P;
    }

    template<typename R>
    polynomial<R> first_chebyshev_polynomial(int n)
    {
        polynomial<R> P1({0,1}),P0({1});
        if(n==0)
            return P0;
        if(n==1)
            return P1;
        else for(int i=2;i<=n;i++)
        {
            std::swap(P0, P1);
            P1 = 2 * P0 * P0 - P1;
        }
        return P1;
    }

    template<typename R>
    std::pair<polynomial<R>,polynomial<R>> fast_chebyshev_polynomial_pair(int n)
    {
        struct chebyshev_poly_pair
        {
            polynomial<R> P0,P1;
            chebyshev_poly_pair& operator*= (const chebyshev_poly_pair &other)
            {
                polynomial<R> A=fast_multiplication(P0,other.P0),B=fast_multiplication(P1,other.P1),C=fast_multiplication(P0,other.P1),D=fast_multiplication(P1,other.P0);
                static polynomial<R> X2({-1,0,1});
                P0=A+B*X2;
                P1=C+D;
                return *this;
            }

            chebyshev_poly_pair operator*(const chebyshev_poly_pair &other) const
            {
                chebyshev_poly_pair P=*this;
                P*=other;
                return P;
            }
            operator std::pair<polynomial<R>,polynomial<R>>()
            {
                return {P0,P1};
            }
        };
        chebyshev_poly_pair P;
        P.P0.data()={0,1};
        P.P1.data()={1};
        return pow(P,n);

    }

    template<typename R>
    polynomial<R> fast_first_chebyshev_polynomial(int n)
    {
        return fast_chebyshev_polynomial_pair<R>(n).first;
    }

    template<typename R>
    polynomial<R> fast_second_chebyshev_polynomial(int n)
    {
        return fast_chebyshev_polynomial_pair<R>(n+1).second;
    }

    template<typename R>
    polynomial<R> cyclotomic_polynomial_prime(int p)
    {
        std::vector<R> P(p-1,1);
        return P;
    }
}


#endif //CPLIBRARY_SPECIAL_POLYNOMIALS_H

//
// Created by ramizouari on 12/01/24.
//

#ifndef CPLIBRARY_GENERATING_H
#define CPLIBRARY_GENERATING_H
#include "signals/fft.h"
#include "signals/ntt.h"
#include "polynomial/fast_polynomial.h"
#include "factorial.h"
#include "polynomial/formal_series.h"
namespace cp
{
    template<cp::integer m>
    struct ogf
    {
        virtual ~ogf()= default;
        std::vector<cp::cyclic<m>> P;
        ogf(std::vector<cp::cyclic<m>> _P): P(_P)
        {
        }

        ogf operator+(const ogf &G) const
        {
            std::vector<cp::cyclic<m>> res;
            for(int i=0;i<P.size();i++)
                res.push_back(P[i]+G.P[i]);
            return ogf(res);
        }
        ogf& operator+=(const ogf &G)
        {
            for(int i=0;i<P.size();i++)
                P[i]+=G.P[i];
            return *this;
        }

        operator const std::vector<cp::cyclic<m>> &() const
        {
            return P;
        }

        operator cp::polynomial<cp::cyclic<m>> () const
        {
            return cp::polynomial<cp::cyclic<m>>(P);
        }

        virtual cp::cyclic<m> weight(cp::integer n) const
        {
            return n<P.size()?P[n]:0;
        }


        const std::vector<cp::cyclic<m>> &data() const
        {
            return P;
        }

        std::vector<cp::cyclic<m>> &data()
        {
            return P;
        }
    };

    template<cp::integer m>
    using ordinary_generating_function=ogf<m>;

    template<integer m>
    ogf<m> multiplication(const ogf<m> &A, const ogf<m> &B)
    {
        return ogf<m>(fast_multiplication(A.data(), B.data()));
    }

    template<cp::integer m>
    ogf<m> power(const ogf<m> &A, integer n)
    {
        if(n==0)
            return 1;
        else
        {
            auto Q=power(A,n/2);
            Q=multiplication(Q,Q);
            if(n%2)
                Q=multiplication(Q,A);
            return Q;
        }
    }

    template<integer m>
    ogf<m> multiplication(const ogf<m> &A, const ogf<m> &B, binary_operation_ptr<std::vector<cyclic<m>>> Op)
    {
        return ogf<m>(Op(A, B));
    }

    template<integer m>
    ogf<m> exp(const ogf<m> &A)
    {
        return cp::formal_exp(cp::polynomial(A.data()),A.data().size()).data();
    }

    template<integer m>
    ogf<m> sequence(const ogf<m> &A)
    {
        int n=A.data().size();
        polynomial<cyclic<m>> P(A);
        P.data().resize(n);
        return formal_inv(cyclic<m>{1}-P,n).data();
    }

    template<integer m>
    ogf<m> multiset(const ogf<m> &A)
    {
        int n=A.data().size();
        polynomial<cyclic<m>> P(A);
        P.data().resize(n);
        for(integer i=2;i<n;i++) for(integer j=1;j<n && i*j < n;j++)
            P[i*j]+=A.data()[j]/i;
        return formal_exp(P,P.data().size()).data();
    }

    template<integer m>
    ogf<m> multiset(const ogf<m> &A, std::vector<integer> &I)
    {
        int n=A.data().size();
        polynomial<cyclic<m>> P(A);
        P.data().resize(n);
        for(integer i=2;i<n;i++) for(integer j=1;j<n && i*j < n;j++)
                P[i*j]+=A.data()[j]*I[i];
        return formal_exp(P,P.data().size()).data();
    }

    template<integer m>
    ogf<m> set(const ogf<m> &A)
    {
        int n=A.data().size();
        polynomial<cyclic<m>> P(A);
        P.data().resize(n);
        bool invert=true;
        for(integer i=2;i<n;i++,invert=!invert) for(integer j=1;j<n && i*j < n;j++)
        {
            if(invert)
                P[i*j]-=A.data()[j]/i;
            else
                P[i*j]+=A.data()[j]/i;
        }
        return formal_exp(P,P.data().size()).data();
    }

    template<integer m>
    ogf<m> set(const ogf<m> &A, std::vector<integer> &I)
    {
        int n=A.data().size();
        polynomial<cyclic<m>> P(A);
        P.data().resize(n);
        bool invert=true;
        for(integer i=2;i<n;i++,invert=!invert) for(integer j=1;j<n && i*j < n;j++)
            {
                if(invert)
                    P[i*j]-=A.data()[j]*I[i];
                else
                    P[i*j]+=A.data()[j]*I[i];
            }
        return formal_exp(P,P.data().size()).data();
    }

}
#endif //CPLIBRARY_GENERATING_H

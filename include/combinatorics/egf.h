//
// Created by ramizouari on 12/01/24.
//

#ifndef CPLIBRARY_EGF_H
#define CPLIBRARY_EGF_H

#include "signals/fft.h"
#include "signals/ntt.h"
#include "polynomial/fast_polynomial.h"
#include "factorial.h"
#include "polynomial/formal_series.h"

namespace cp
{
    template<cp::integer m>
    struct egf
    {
        std::vector<cp::cyclic<m>> P;
        std::shared_ptr<abstract_factorial<cp::cyclic<m>>> factorial;
        std::shared_ptr<abstract_inverse_factorial<cp::cyclic<m>>> inverse_factorial;
        egf(const std::vector<cp::cyclic<m>> &_P, std::shared_ptr<abstract_factorial<cp::cyclic<m>>> _factorial, bool normalized=true,
            std::shared_ptr<abstract_inverse_factorial<cp::cyclic<m>>> _inverse_factorial=nullptr): P(_P)
        {
            factorial=_factorial;
            inverse_factorial=_inverse_factorial;
            if(!normalized)
            {
                if(inverse_factorial) for(int i=0;i<P.size();i++)
                        P[i]*=inverse_factorial->inverse_factorial(i);
                else for(int i=0;i<P.size();i++)
                        P[i]/=factorial->factorial(i);
            }
        }

        cp::cyclic<m> weight(cp::integer n) const
        {
            return n<P.size()?P[n]*factorial->factorial(n):0;
        }

        egf operator+(const egf &G) const
        {
            std::vector<cp::cyclic<m>> res;
            for(int i=0;i<P.size();i++)
                res.push_back(P[i]+G.P[i]);
            return ogf(res);
        }
        egf& operator+=(const egf &G)
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
    using exponential_generating_function=egf<m>;

    template<integer m>
    egf<m> exp_series(integer x,integer n,std::shared_ptr<abstract_factorial<cp::cyclic<m>>> F)
    {
        std::vector<cyclic<m>> P(n);
        P[0]=1;
        for(int i=1;i<n;i++)
            P[i]=P[i-1]*x/F->factorial(i);
        return egf<m>(P,F);
    }

    template<integer m>
    egf<m> exp_series(integer x,integer n,std::shared_ptr<abstract_factorial<cp::cyclic<m>>> F,std::shared_ptr<abstract_inverse_factorial<cp::cyclic<m>>> F_inv)
    {
        std::vector<cyclic<m>> P(n);
        P[0]=1;
        for(int i=1;i<n;i++)
            P[i] = P[i - 1] * x;
        for(int i=0;i<n;i++)
            P[i]*=F_inv->inverse_factorial(i);
        return egf<m>(P,F,true,F_inv);
    }



    template<integer m>
    egf<m> multiplication(const egf<m> &A, const egf<m> &B)
    {
        auto Z=fast_multiplication(A.data(), B.data());
        Z.resize(std::max(A.P.size(),B.data().size()));
        return egf<m>(Z,A.factorial,true,A.inverse_factorial);
    }

    template<integer m>
    egf<m> multiplication(const egf<m> &A, const egf<m> &B, binary_operation_ptr<std::vector<cyclic<m>>> Op)
    {
        return egf<m>(Op(A, B));
    }

    template<cp::integer m>
    egf<m> power(const egf<m> &A, integer n)
    {
        if(n==0)
            return egf<m>({1},A.factorial,true,A.inverse_factorial);
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
    egf<m> exp(const egf<m> &A)
    {
        return egf<m>(formal_exp(cp::polynomial(A.data()),A.data().size()).data(),A.factorial,true,A.inverse_factorial);
    }

    template<integer m>
    egf<m> sequence(const egf<m> &A)
    {
        int n=A.data().size();
        polynomial<cyclic<m>> P(A);
        P.data().resize(n);
        return egf<m>(formal_inv(cyclic<m>{1}-P,n).data(),A.factorial,true,A.inverse_factorial);
    }


    template<integer m>
    egf<m> set(const egf<m> &A)
    {
        return egf<m>(formal_exp(polynomial(A.data()),A.data().size()).data(),A.factorial,true,A.inverse_factorial);
    }

}

#endif //CPLIBRARY_EGF_H

//
// Created by ramizouari on 11/01/24.
//

#ifndef CPLIBRARY_FACTORIAL_H
#define CPLIBRARY_FACTORIAL_H
#include <concepts>
#include "algebra/abstract_algebra.h"
#include "nt/modular/fixed.h"
#include "nt/modular_functions.h"

namespace cp
{
    template<typename T>
    struct abstract_factorial
    {
        virtual ~abstract_factorial()= default;
        virtual T factorial(integer n) const=0;
        T operator()(integer n) const
        {
            return factorial(n);
        }
    };

    template<typename T>
    struct abstract_inverse_factorial
    {
        virtual T inverse_factorial(integer n) const=0;
        T operator()(integer n) const
        {
            return inverse_factorial(n);
        }
    };

    template<typename T>
    struct naive_factorial : public abstract_factorial<T>
    {
        T factorial(integer n) const override
        {
            T res=1;
            for(int i=2;i<=n;i++)
                res*=i;
            return res;
        }
    };

    template<typename T>
    struct precomputed_factorial : public abstract_factorial<T>
    {
        std::vector<T> F;
        explicit precomputed_factorial(T n)
        {
            F.resize(n+1);
            F[0]=1;
            for(T i=1;i<=n;i++)
                F[i]=F[i-1]*i;
        }
        T factorial(integer n)  const override
        {
            return F[n];
        }
    };

    template<integer n>
    struct precomputed_factorial<cp::cyclic<n>> : public abstract_factorial<cp::cyclic<n>>,public abstract_inverse_factorial<cp::cyclic<n>>
    {
        std::vector<cp::cyclic<n>> F,F_inv;

        explicit precomputed_factorial(integer m, const std::vector<integer> &I):F(m+1),F_inv(m+1)
        {
            F[0]=1;
            F_inv[0]=1;
            for(integer i=1;i<=m;i++)
            {
                F[i] = F[i - 1] * i;
                F_inv[i] = F_inv[i - 1] * I[i];
            }
        }

        explicit precomputed_factorial(integer m): precomputed_factorial(m,cp::inverse_table(m,cp::cyclic<n>::modulus()))
        {

        }

        cp::cyclic<n> factorial(integer m)  const override
        {
            return F[m];
        }

        cp::cyclic<n> inverse_factorial(integer m) const override
        {
            return F_inv[m];
        }
        using abstract_factorial<cp::cyclic<n>>::operator();
    };

    template<>
    struct precomputed_factorial<cp::cyclic<cp::dynamic_modulus>> : public abstract_factorial<cp::cyclic<cp::dynamic_modulus>>,public abstract_inverse_factorial<cp::cyclic<cp::dynamic_modulus>>
    {
        std::vector<cp::cyclic<cp::dynamic_modulus>> F,F_inv;
        explicit precomputed_factorial(integer m,integer p, const std::vector<integer> &I):F(m+1),F_inv(m+1)
        {
            F[0]={1,p};
            F_inv[0]={1,p};
            for(integer i=1;i<=m;i++)
            {
                F[i] = F[i - 1] * i;
                F_inv[i] = F_inv[i - 1] * I[i];
            }
        }

        explicit precomputed_factorial(integer m,integer p): precomputed_factorial(m,p,cp::inverse_table(m,p))
        {

        }

        cp::cyclic<cp::dynamic_modulus> factorial(integer m)  const override
        {
            return F[m];
        }

        cp::cyclic<cp::dynamic_modulus> inverse_factorial(integer m) const override
        {
            return F_inv[m];
        }
        using abstract_factorial::operator();
    };
}

#endif //CPLIBRARY_FACTORIAL_H

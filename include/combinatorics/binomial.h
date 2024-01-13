//
// Created by ramizouari on 12/01/24.
//

#ifndef CPLIBRARY_BINOMIAL_NUMBERS_H
#define CPLIBRARY_BINOMIAL_NUMBERS_H
#include "factorial.h"
namespace cp
{
    template<typename T>
    struct abstract_nCr
    {
        virtual T nCr(integer n,integer r) const=0;
        T operator()(integer n,integer r) const
        {
            return nCr(n,r);
        }
        virtual ~abstract_nCr()= default;
    };

    template<typename T>
    struct naive_nCr : public abstract_nCr<T>
    {
        T nCr(integer n,integer r) const override
        {
            if(r>n)
                return 0;
            T res=1;
            for(int i=0;i<r;i++)
                res*=(n-i);
            for(int i=1;i<=r;i++)
                res/=i;
            return res;
        }
    };

    template<typename T>
    struct nCr_t : public abstract_nCr<T>
    {
        std::shared_ptr<abstract_factorial<T>> factorial;
        std::shared_ptr<abstract_inverse_factorial<T>> inverse_factorial;
        explicit nCr_t(std::shared_ptr<abstract_factorial<T>> _factorial,std::shared_ptr<abstract_inverse_factorial<T>> _inverse_factorial=nullptr)
        {
            factorial=_factorial;
            inverse_factorial=_inverse_factorial;
        }
        T nCr(integer n,integer r) const override
        {
            if(r>n)
                return 0;
            if(inverse_factorial)
                return (*factorial)(n)*(*inverse_factorial)(r)*(*inverse_factorial)(n-r);
            else
                return (*factorial)(n)/(*factorial)(r)/(*factorial)(n-r);
        }
    };

    template<typename T>
    struct precomputed_nCr : public abstract_nCr<T>
    {
        std::vector<std::vector<T>> C;
        explicit precomputed_nCr(T n)
        {
            C.resize(n+1);
            C[0][0]=1;
            for(int i=0;i<=n;i++) for(int j=0;j<=i;j++)
                {
                    [[likely]]
                    if(i>0)
                        C[i][j]+=C[i-1][j];
                    [[likely]]
                    if(i>0 && j>0)
                        C[i][j]+=C[i-1][j-1];
                }
        }
        T nCr(integer n,integer r) const override
        {
            return C[n][r];
        }
    };
}
#endif //CPLIBRARY_BINOMIAL_H

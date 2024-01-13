//
// Created by ramizouari on 12/01/24.
//

#ifndef CPLIBRARY_PARTITIONS_H
#define CPLIBRARY_PARTITIONS_H

#include "algebra/abstract_algebra.h"

namespace cp
{
    template<typename T>
    struct abstract_partition_number
    {
        virtual ~abstract_partition_number()= default;
        virtual T partition_number(integer n) const=0;
        T operator()(integer n) const
        {
            return partition_number(n);
        }
    };

    template<typename T>
    struct precomputed_partition_number : public abstract_partition_number<T>
    {
        std::vector<T> P;
        explicit precomputed_partition_number(integer n)
        {
            P.resize(n+1);
            P[0]=1;
            for(int i=1;i<=n;i++) for(int j=1;j<=i;j++)
                    P[i]+=P[i-j];
        }
        T partition_number(integer n) const override
        {
            return P[n];
        }
    };

    template<typename T>
    struct euler_partition_number : public abstract_partition_number<T>
    {
        std::vector<T> P;
        explicit euler_partition_number(integer n)
        {
            P.resize(n+1);
            P[0]=1;
            for(int i=1;i<=n;i++)
            {
                T res=0;
                for(integer k=1;;k++)
                {
                    integer gk=k*(3*k-1)/2;
                    if(gk>i) break;
                    res+=P[i-gk];
                }
                for(integer k=1;;k++)
                {
                    integer gk=k*(3*k+1)/2;
                    if(gk>i) break;
                    res+=P[i-gk];
                }
                P[i]=res;
            }
        }
        T partition_number(integer n) const override
        {
            return P[n];
        }
    };
}

#endif //CPLIBRARY_PARTITIONS_H

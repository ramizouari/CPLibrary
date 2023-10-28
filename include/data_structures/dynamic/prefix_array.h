//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_DYN_PREFIX_ARRAY_H
#define CPLIBRARY_DYN_PREFIX_ARRAY_H
#include <vector>
#include "algebra/binary_operation.h"
#include <memory>

namespace data_structures::dynamic
{

    template<typename R>
    struct prefix_array
    {
        std::vector<R> A;
        std::vector<R> P;
        invertible_binary_operation_ptr<R> F;
        prefix_array(const std::vector<R> &_A, std::shared_ptr<binary_operation<R>> _F, std::shared_ptr<invertible_operation<R>> _I):A(_A),P(_A.size()+1),F(_F,_I)
        {
            P[0]=F->neutral_element();
            for(int i=0;i<A.size();i++)
                P[i+1]=F(P[i],A[i]);
        }

        R query(int l,int r)
        {
            return F(F.inv(P[l]),P[r]);
        }

        void update(int i,R u)
        {
            A[i]=u;
            for(int j=i+1;j<P.size();j++)
                P[j]=F(P[j-1],A[j-1]);
        }
    };
}

#endif //CPLIBRARY_PREFIX_ARRAY_H

//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_DYN_SPARSE_ARRAY_H
#define CPLIBRARY_DYN_SPARSE_ARRAY_H
#include <vector>
#include "algebra/binary_operation.h"
#include <memory>
#include "algebra/bits.h"

namespace cp::data_structures::dynamic
{
    template<typename T>
    struct sparse_array
    {
        int n,h;
        std::vector<std::vector<T>> S;
        binary_operation_ptr<T> F;
    public:
        sparse_array(const std::vector<T>&A, std::shared_ptr<binary_operation<T>> _F):n(bit_ceil(A.size())),h(bit_log(n)),S(h+1),F(_F)
        {
            int r=1;
            for(int i=h;i>=0;i--,r*=2)
                S[i].resize(n-r+1,F->neutral_element());
            for(int i=0;i<A.size();i++)
                S[h][i]=A[i];
            r=1;
            for(int i=h-1;i>=0;i--,r*=2) for(int j=0;j<=n-2*r;j++)
                    S[i][j]=F(S[i+1][j],S[i+1][j+r]);
        }

        T query(int l,int r) const
        {
            if(l>=r)
                return F->neutral_element();
            auto d=r-l;
            auto s=bit_floor(d);
            auto b=bit_log(s);
            return F(S[h-b][l],S[h-b][r-s]);
        }
    };
}

#endif //CPLIBRARY_SPARSE_ARRAY_H

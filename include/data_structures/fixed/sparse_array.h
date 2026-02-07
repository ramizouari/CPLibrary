//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_FIXED_SPARSE_ARRAY_H
#define CPLIBRARY_FIXED_SPARSE_ARRAY_H
#include <vector>
#include "algebra/bits.h"
namespace cp::data_structures::fixed
{
    template<typename O>
    struct sparse_array
    {
        using T= O::type;
        using type = T;
        using value_type = T;
        using key_type = int;
        using binary_operation = O;
        inline static O F=O();
        int n,h;
        std::vector<std::vector<T>> S;

        sparse_array(const std::vector<T>&A):n(bit_ceil(A.size())),h(bit_log(n)),S(h+1)
        {
            int r=1;
            for(int i=h;i>=0;i--,r*=2)
                S[i].resize(n-r+1,O::neutral);
            for(int i=0;i<A.size();i++)
                S[h][i]=A[i];
            r=1;
            for(int i=h-1;i>=0;i--,r*=2) for(int j=0;j<=n-2*r;j++)
                    S[i][j]=F(S[i+1][j],S[i+1][j+r]);
        }

        T query(int l,int r) const
        {
            if(l>=r)
                return O::neutral;
            auto d=r-l;
            auto s=bit_floor(d);
            auto b=bit_log(s);
            return F(S[h-b][l],S[h-b][r-s]);
        }

        std::vector<T> data() const
        {
            return S[h];
        }

        std::span<const T> span() const
        {
            return S[h];
        }
    };


}
#endif //CPLIBRARY_SPARSE_ARRAY_H

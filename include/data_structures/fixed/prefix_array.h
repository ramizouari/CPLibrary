//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_FIXED_PREFIX_ARRAY_H
#define CPLIBRARY_FIXED_PREFIX_ARRAY_H
#include <vector>
namespace cp::data_structures::fixed
{
    template<typename O>
    struct prefix_array
    {
        using R=typename O::type;
        using type=typename O::type;
        std::vector<R> A;
        std::vector<R> P;
        inline static O F=O();
        prefix_array(const std::vector<R> &_A):A(_A),P(_A.size()+1)
        {
            P[0]=O::neutral;
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

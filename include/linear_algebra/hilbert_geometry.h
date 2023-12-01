//
// Created by ramizouari on 09/12/22.
//

#ifndef CPLIBRARY_HILBERT_GEOMETRY_H
#define CPLIBRARY_HILBERT_GEOMETRY_H


#include "decomposition.h"

namespace cp::linalg
{
    template<typename IK>
    IK inner_product(const d_vector<IK>& u,const d_vector<IK> &v)
    {
        IK r=0;
        int n=std::min(u.dim(),v.dim());
        for(int i=0;i<n;i++)
            r+=conj<IK,IK>(u)*v;
        return r;
    }

    template<typename IK>
    void reflect_inplace(d_vector<IK> &x,const d_vector<IK> &u)
    {
        x-=2* unit_proj(u,x);
    }

    template<typename IK>
    d_vector<IK> reflect(const d_vector<IK> &x,const d_vector<IK> &u)
    {
        auto y=x;
        reflect_inplace(y,u);
        return y;
    }

    template<typename IK>
    void rotate_inplace(d_vector<IK>& x,const d_vector<IK> &u,const d_vector<IK>&v)
    {

    }
}


#endif //CPLIBRARY_HILBERT_GEOMETRY_H

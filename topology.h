#ifndef __TOPOLOGY_H__
#define __TOPOLOGY_H__
#include "abstract_algebra.h"
#include "linear_algebra.h"

template<typename E>
struct metric_t
{
    virtual real metric(const E&u,const E&v)const =0;
};

template<typename E>
struct norm_t :public metric_t<E>
{
    virtual real norm(const E&u)const =0;
    real metric(const E&u,const E&v)const  override
    {
        return norm(v-u);
    }
};

template<typename E>
struct inner_product_t:public norm_t<E>
{
    real norm(const E&u) const override
    {
        return std::sqrt(inner_product(u,u));
    }

    virtual real inner_product(const E&u,const E&v)const  =0;
};

template<typename E>
struct L2_inner_product :public inner_product_t<E>
{
    real inner_product(const E&u,const E&v) const
    {
        auto m=std::min(u.dim(),v.dim());
        real R=0;
        for(int i=0;i<m;i++)
            R+=u[i]*v[i];
        return R;
    }
};

template<typename E>
struct L1_norm :public norm_t<E>
{
    real norm(const E&u) const
    {
        real R=0;
        for(int i=0;i<u.dim();i++)
            R+=std::abs(u[i]);
        return R;
    }
};

template<typename E>
struct L_inf_norm :public norm_t<E>
{
    real norm(const E&u) const
    {
        real R=0;
        for(int i=0;i<u.dim();i++)
            R=std::max(R,std::abs(u[i]));
        return R;
    }
};

struct real_inner_product :public inner_product_t<real>
{
    real inner_product(const real &u,const real &v) const override
    {
        return u*v;
    }
};

template<>
struct L2_inner_product<real>:public real_inner_product{};
template<>
struct L_inf_norm<real>:public real_inner_product{};
template<>
struct L1_norm<real>:public real_inner_product{};
#endif
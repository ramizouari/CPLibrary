//
// Created by ASUS on 01/12/2021.
//
#ifndef __OPERATION_H__
#define __OPERATION_H__
#include <numeric>
#include "abstract_algebra.h"

template<typename T>
struct binary_operation
{
    virtual T operator()(const T&a,const T&b) const=0;
    template<typename H0,typename ...H>
    T operator()(const H0&a,const H&... b)
    {
        return this->operator()(a,this->operator()(b...));
    }
};

template<typename T>
struct plus_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return a+b;
    }

    T inv(const T&a) const
    {
        return -a;
    }
    inline static T neutral=T();
};

template<typename T>
struct multiplies_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return a*b;
    }

    inline static T neutral=T(1);
};

template<typename T>
struct field_multiplies_t:public multiplies_t<T>
{
    T inv(const T&a) const
    {
        return a.inv();
    }
};

template<>
struct field_multiplies_t<real>:public multiplies_t<real>
{
    real inv(const real& a)const
    {
        return 1./a;
    }
};

template<>
struct field_multiplies_t<IC>:public multiplies_t<IC>
{
    IC inv(const IC& a)const
    {
        return 1./a;
    }
};

template<typename T>
struct max_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return std::max(a,b);
    }

    inline static T neutral=0;
};

template<typename T>
struct min_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return std::min(a,b);
    }

    inline static T neutral=1e9;
};

template<typename T>
struct gcd_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return std::gcd(a,b);
    }

    inline static T neutral=0;
};

template<typename T>
struct lcm_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return std::lcm(a,b);
    }

    inline static T neutral=1;
};

template<typename T>
struct xor_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return a^b;
    }

    T inv(const T&a) const
    {
        return a;
    }

    inline static T neutral=0;
};

template<typename T>
struct and_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return a&b;
    }

    inline static T neutral=static_cast<T>(-1);
};

template<typename T>
struct or_t:public binary_operation<T>
{
    T operator()(const T&a,const T&b) const
    {
        return a|b;
    }

    inline static T neutral=0;
};

#endif
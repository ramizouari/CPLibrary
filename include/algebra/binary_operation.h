//
// Created by ASUS on 01/12/2021.
//
#ifndef __OPERATION_H__
#define __OPERATION_H__
#include <numeric>
#include "abstract_algebra.h"
#include <memory>

template<typename T>
struct binary_operation
{
    using type=T;
    template<typename H0,typename ...H>
    T operator()(const H0&a,const H&... b) const
    {
        if constexpr (sizeof...(b) == 0)
            return a;
        else return reduce(a,this->operator()(b...));
    }
    virtual T reduce(const T& a, const T& b) const = 0;
    virtual T neutral_element() const
    {
        return T{};
    }
};

template<typename T>
struct invertible_operation
{
    virtual T inv(const T& a) const = 0;
};

template<typename T>
struct monoid_plus_t:public binary_operation<T> {
    T reduce(const T &a, const T &b) const override {
        return a + b;
    }
    inline static T neutral{};
};

template<typename T>
struct plus_t:public monoid_plus_t<T>,public invertible_operation<T>
{
    T inv(const T&a) const override
    {
        return -a;
    }
};

template<typename T>
struct multiplies_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return a*b;
    }

    inline static T neutral=T(1);
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct field_multiplies_t:public multiplies_t<T>,public invertible_operation<T>
{
    T inv(const T&a) const
    {
        return a.inv();
    }

};

template<>
struct field_multiplies_t<real>:public multiplies_t<real>,public invertible_operation<real>
{
    real inv(const real& a)const
    {
        return 1./a;
    }
};

template<>
struct field_multiplies_t<IC>:public multiplies_t<IC>,public invertible_operation<IC>
{
    IC inv(const IC& a)const
    {
        return IC(1)/a;
    }
};

template<typename T>
struct max_t:public binary_operation<T>
{
    T e;
    explicit max_t(T _e):e(_e){}
    max_t(): max_t(T{}){}
    T reduce(const T&a,const T&b) const override
    {
        return std::max(a,b);
    }

    inline static T neutral{0};
    T neutral_element() const override
    {
        return e;
    }
};

template<typename T>
struct min_t:public binary_operation<T>
{
    T e;
    explicit min_t(T _e):e(_e){}
    min_t(): min_t(T{}){}

    T reduce(const T&a,const T&b) const override
    {
        return std::min(a,b);
    }

    inline static T neutral{};

    T neutral_element() const override
    {
        return e;
    }
};

template<typename T>
struct gcd_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return gcd(a,b);
    }

    inline static T neutral{0};
};

template<typename T>
struct lcm_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return lcm(a,b);
    }

    inline static T neutral{1};
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct xor_t:public binary_operation<T>,public invertible_operation<T>
{
    T reduce(const T&a,const T&b) const
    {
        return a^b;
    }

    T inv(const T&a) const
    {
        return a;
    }

    inline static T neutral{};
};

template<typename T>
struct and_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return a&b;
    }

    inline static T neutral=static_cast<T>(-1);
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct or_t:public binary_operation<T>
{
    T reduce(const T&a,const T&b) const override
    {
        return a|b;
    }

    inline static T neutral{};
};

template<typename T>
struct logical_and_t :public binary_operation<T>
{
    T reduce(const T& a, const T& b) const override
    {
        return a && b;
    }

    inline static T neutral{true};
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct logical_or_t :public binary_operation<T>
{
    T reduce(const T& a, const T& b) const override
    {
        return a || b;
    }

    inline static T neutral{false};
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
struct logical_xor_t :public binary_operation<T>,public invertible_operation<T>
{
    T reduce(const T& a, const T& b) const override
    {
        return !a && b || a && !b;
    }
    T inv(const T&a) const
    {
        return !a;
    }
    inline static T neutral{false};
    T neutral_element() const override
    {
        return neutral;
    }
};

template<typename T>
class binary_operation_ptr
{
    std::shared_ptr<binary_operation<T>> op;
public:
    binary_operation_ptr(std::shared_ptr<binary_operation<T>> value): op(value){}
    template<typename ...H>
    auto operator()(const H&... h) const
    {
        return op->operator()(h...);
    }

    auto neutral_element() const
    {
        return op->neutral_element();
    }
};

template<typename T>
class invertible_binary_operation_ptr : public binary_operation_ptr<T>
{
    std::shared_ptr<invertible_operation<T>> inverter;
public:
    invertible_binary_operation_ptr(std::shared_ptr<binary_operation<T>> b,
                                    std::shared_ptr<invertible_operation<T>> I): binary_operation_ptr<T>(b),inverter(I){}
    using binary_operation_ptr<T>::operator();
    using binary_operation_ptr<T>::neutral_element;
    auto inv(const T& a) const
    {
        return inverter->inv(a);
    }
};

#endif
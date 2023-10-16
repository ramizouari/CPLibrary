//
// Created by ramizouari on 09/12/22.
//

#ifndef CPLIBRARY_VECTOR_H
#define CPLIBRARY_VECTOR_H

#include <vector>
#include <cstddef>
#include <array>

struct v_shape
{
    int n;
};

/**
 * @brief Dynamic Vector
 * @detail Dynamic Vector is a vector in the mathematical sense,
 * @formal It is the union of R^k for all k, where k is the dimension of the vector.
 * <ul>
 * <li> Addition between 2 vectors are defined with respect to the first vector's shape <br>
 * <li> for all k, the set of all vectors of shape k is a vector space
 * </ul>
 * @Requirements
 * R is a commutative ring
 * */
template<typename R>
class d_vector
{
    std::vector<R> u;
public:
    using base_field=R;
    using base_ring=R;
    inline static int n=0;
    d_vector():u(n){}
    d_vector(std::vector<R> _u):u(std::move(_u)){}
    d_vector(v_shape shape):u(shape.n){}

    bool operator==(const d_vector<R>& other) const
    {
        return u==other.u;
    }
    auto dim() const
    {
        return u.size();
    }

    auto& operator[](int k)
    {
        return u[k];
    }

    const auto& operator[](int k) const
    {
        return u[k];
    }

    auto& operator+=(const d_vector &o)
    {
        for(int i=0;i<dim();i++)
            u[i]+=o.u[i];
        return *this;
    }

    auto& operator-=(const d_vector &o)
    {
        for(int i=0;i<dim();i++)
            u[i]-=o.u[i];
        return *this;
    }

    auto& operator*=(R k)
    {
        for(auto &s:u)
            s*=k;
        return *this;
    }

    auto operator+(const d_vector &o) const
    {
        auto v=*this;
        return v+=o;
    }

    auto operator-(const d_vector &o) const
    {
        auto v=*this;
        return v-=o;
    }

    auto operator-() const
    {
        auto v=*this;
        for(auto &s:v.u)
            s=-s;
        return v;
    }

    auto& operator/=(R k)
    {
        for(auto &s:u)
            s/=k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto v=*this;
        return v/=k;
    }

    auto begin() {
        return u.begin();
    }

    auto begin() const
    {
        return u.cbegin();
    }

    auto end()
    {
        return u.end();
    }

    auto end() const
    {
        return u.cend();
    }
};

template<typename R>
auto operator*(const R&k,const d_vector<R>& u)
{
    auto v=u;
    return v*=k;
}

/**
 * @brief Static Vector:
 * @tparam R is the base field
 * @tparam n is the dimension of the vector space
 * @details It is a member of an R-vector space E where dim(E)= n
 * @Requirements
 * <strong>R</strong> is a commutative ring. <br>
 * @Formal <strong>E</strong> is an <strong>R</strong>-module, and it is a vector space only if <strong>R</strong> is a field. <br>
 * In fact, the name s_vector is used for consistency with the computer science's name.
 */

template<typename R,int n>
class s_vector
{
    std::array<R,n> u;
public:
    using base_field=R;
    using base_ring=R;
    inline static constexpr int dim()
    {
        return n;
    }

    s_vector()
    {
        for(int i=0;i<n;i++)
            u[i]=0;
    }

    s_vector(std::array<R,n>_u):u(std::move(_u)){}

    bool operator==(const s_vector&) const = default;

    auto& operator[](int k)
    {
        return u[k];
    }

    const auto& operator[](int k) const
    {
        return u[k];
    }

    auto& operator+=(const s_vector &o)
    {
        auto r=std::min(dim(),o.dim());
        for(int i=0;i<r;i++)
            u[i]+=o.u[i];
        return *this;
    }

    auto& operator-=(const s_vector &o)
    {
        auto r=std::min(dim(),o.dim());
        for(int i=0;i<r;i++)
            u[i]-=o.u[i];
        return *this;
    }

    auto& operator*=(R k)
    {
        for(auto &s:u)
            s*=k;
        return *this;
    }

    auto operator+(const s_vector &o) const
    {
        auto v=*this;
        return v+=o;
    }

    auto operator-(const s_vector &o) const
    {
        auto v=*this;
        return v-=o;
    }

    auto operator-() const
    {
        auto v=*this;
        for(auto &s:v.u)
            s=-s;
        return v;
    }

    auto& operator/=(R k)
    {
        for(auto &s:u)
            s/=k;
        return *this;
    }

    auto operator/(R k) const
    {
        auto v=*this;
        return v/=k;
    }

    auto begin()
    {
        return u.begin();
    }

    auto begin() const
    {
        return u.cbegin();
    }

    auto end()
    {
        return u.end();
    }

    auto end() const
    {
        return u.cend();
    }

    template <size_t k>
    auto& get()& {
        return u[k];
    }

    template <size_t k>
    const auto& get() const& {
        return u[k];
    }

    template <size_t k>
    auto&& get() const&& {
        return u[k];
    }

    template <size_t k>
    auto&& get() && {
        return u[k];
    }
};

template<typename R,int n>
auto operator*(const R&k,const s_vector<R,n>& u)
{
    auto v=u;
    return v*=k;
}

namespace std
{
    template<typename R,int n>
    struct tuple_size<s_vector<R, n>> : std::integral_constant<size_t, n>{};
    template<size_t k,typename R,int n>
    struct tuple_element<k, s_vector<R, n>>
    {
        using type = R;
    };
}

#endif //CPLIBRARY_VECTOR_H

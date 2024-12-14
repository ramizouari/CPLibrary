//
// Created by ramizouari on 09/12/22.
//

#ifndef CPLIBRARY_VECTOR_H
#define CPLIBRARY_VECTOR_H

#include <vector>
#include <cstddef>
#include <array>
#include "utils.h"
#include "view.h"
#include "algebra/structures.h"

namespace cp::linalg
{
    using v_shape = std::array<std::size_t, 1>;
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
    template<ring R>
    struct vector : tensor_view<R,1>
    {
        std::vector<R> u;
        using base_field=R;
        using base_ring=R;
        vector() = default;
        vector(std::vector<R>&& _u):u(std::move(_u)){}
        vector(const std::vector<R>& _u):u(_u){}
        vector(std::size_t n,size_tag_t):u(n){}

        bool operator==(const vector& other) const = default;

        std::size_t dim() const
        {
            return u.size();
        }

        std::array<std::size_t,1> shape() const override {
            return {u.size()};
        }

        std::size_t size() const override
        {
            return u.size();
        }

        R& at(std::array<std::size_t,1> i) override
        {
            return u[i.front()];
        }

        const R& at(std::array<std::size_t,1> i) const override
        {
            return u[i.front()];
        }

        R& operator[](std::size_t k)
        {
            return u[k];
        }

        const R& operator[](std::size_t k) const
        {
            return u[k];
        }

        vector& operator+=(const vector &O)
        {
            for(int i=0;i<std::min(dim(),O.dim());i++)
                u[i]+=O.u[i];
            return *this;
        }

        vector& operator-=(const vector &O)
        {
            for(int i=0;i<std::min(dim(),O.dim());i++)
                u[i]-=O.u[i];
            return *this;
        }

        vector& operator*=(const R& k)
        {
            for(auto &s:u) s*=k;
            return *this;
        }


        vector operator-() const
        {
            auto v=*this;
            for(auto &s:v.u) s=-s;
            return v;
        }

        vector& operator/=(const R& k)
        {
            for(auto &s:u) s/=k;
            return *this;
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

        auto& data()
        {
            return u;
        }

        const auto& data() const
        {
            return u;
        }
    };

    template<typename Vec,typename R>
    concept ToVector=std::convertible_to<Vec,vector<R>>;

    template<typename Vec,typename R>
    concept ToVectorProper=std::convertible_to<Vec,vector<R>> && !std::same_as<Vec,vector<R>>;


    template<ring R ,ToVector<R> O>
    vector<R> operator+(const vector<R> &A,const O &B)
    {
        auto C=A;
        return C+=B;
    }

    template<ring R ,ToVectorProper<R> O>
    vector<R> operator+(const O & A,const vector<R> & B)
    {
        vector<R> C=A;
        return C+=B;
    }

    template<ring R ,ToVector<R> O>
    vector<R> operator-(const vector<R> &A,const O &B)
    {
        auto C=A;
        return C-=B;
    }

    template<ring R ,ToVectorProper<R> O>
    vector<R> operator-(const O & A,const vector<R> & B)
    {
        vector<R> C=A;
        return C-=B;
    }

    template<ring R ,std::convertible_to<R> O>
    vector<R> operator*(const O &k,const vector<R> &A)
    {
        auto C=A;
        return C*=k;
    }

    template<ring R ,std::convertible_to<R> O>
    vector<R> operator*(const vector<R> &A,const O &k)
    {
        auto C=A;
        return C*=k;
    }

    template<ring R ,std::convertible_to<R> O>
    vector<R> operator/(const vector<R> &A,const O &k)
    {
        auto C=A;
        return C/=k;
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

        s_vector(const std::array<R,n> &_u):u(_u){}
        s_vector(std::array<R,n> &&_u):u(std::move(_u)){}

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
}

namespace std
{
    template<typename R,int n>
    struct tuple_size<cp::linalg::s_vector<R, n>> : std::integral_constant<size_t, n>{};
    template<size_t k,typename R,int n>
    struct tuple_element<k, cp::linalg::s_vector<R, n>>
    {
        using type = R;
    };
}

#endif //CPLIBRARY_VECTOR_H

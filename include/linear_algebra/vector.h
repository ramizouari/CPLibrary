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
    template<ring R,std::size_t ext = dynamic_extent>
    struct vector : tensor_view<R,1>
    {
        using index_array=typename tensor_view<R,1>::index_array;
        using container = std::conditional_t<ext==dynamic_extent,std::vector<R>,std::array<R,ext>>;
        container u{};
        using base_field=R;
        using base_ring=R;
        vector() = default;
        vector(container&& _u):u(std::move(_u)){}
        vector(const container& _u):u(_u){}
        vector(std::size_t n,size_tag_t) requires(ext==dynamic_extent) :u(n){}
        vector(std::size_t n,size_tag_t) requires(ext!=dynamic_extent)
        {
            [[unlikely]]
            if (ext!=n) throw std::runtime_error("Vector size not consistent with its template");
        }

        bool operator==(const vector& other) const = default;

        std::size_t dim() const
        {
            return u.size();
        }

        index_array shape() const override {
            return {u.size()};
        }

        std::size_t size() const override
        {
            return u.size();
        }

        R& at(index_array i) override
        {
            return u[i.front()];
        }

        const R& at(index_array i) const override
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

    template<typename Vec,typename R,std::size_t ext = dynamic_extent>
    concept ToVector=std::convertible_to<Vec,vector<R,ext>>;

    template<typename Vec,typename R,std::size_t ext = dynamic_extent>
    concept ToVectorProper=std::convertible_to<Vec,vector<R,ext>> && !std::same_as<Vec,vector<R,ext>>;


    template<ring R,std::size_t ext ,ToVector<R,ext> O>
    vector<R,ext> operator+(const vector<R,ext> &A,const O &B)
    {
        auto C=A;
        return C+=B;
    }

    template<ring R,std::size_t ext ,ToVectorProper<R,ext> O>
    vector<R,ext> operator+(const O & A,const vector<R,ext> & B)
    {
        vector<R,ext> C=A;
        return C+=B;
    }

    template<ring R,std::size_t ext ,ToVector<R,ext> O>
    vector<R,ext> operator-(const vector<R,ext> &A,const O &B)
    {
        auto C=A;
        return C-=B;
    }

    template<ring R,std::size_t ext ,ToVectorProper<R,ext> O>
    vector<R,ext> operator-(const O & A,const vector<R,ext> & B)
    {
        vector<R> C=A;
        return C-=B;
    }

    template<ring R,std::size_t ext,std::convertible_to<R> O>
    vector<R,ext> operator*(const O &k,const vector<R,ext> &A)
    {
        auto C=A;
        return C*=k;
    }

    template<ring R,std::size_t ext ,std::convertible_to<R> O>
    vector<R,ext> operator*(const vector<R,ext> &A,const O &k)
    {
        auto C=A;
        return C*=k;
    }

    template<ring R,std::size_t ext ,std::convertible_to<R> O>
    vector<R,ext> operator/(const vector<R,ext> &A,const O &k)
    {
        auto C=A;
        return C/=k;
    }
}


#endif //CPLIBRARY_VECTOR_H

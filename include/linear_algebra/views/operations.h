//
// Created by ramizouari on 14/12/24.
//

#ifndef CPLIBRARY_VIEW_OPERATIONS_H
#define CPLIBRARY_VIEW_OPERATIONS_H
#include "../view.h"
#include "algebra/structures.h"
#include "functional/zip.h"
#include "linear_algebra/tensor.h"

namespace cp::linalg {
    template<additive_monoid M,std::size_t Rank>
    void add(tensor_view<M,Rank> &c,const tensor_view<M,Rank> &a,const tensor_view<M,Rank> &b) {
        for (auto && [z,x,y] : zip(c,a,b))
            z = x + y;
    }

    template<additive_group M,std::size_t Rank>
    void subtract(tensor_view<M,Rank> &c,const tensor_view<M,Rank> &a,const tensor_view<M,Rank> &b) {
        for (auto && [z,x,y] : zip(c,a,b))
            z = x - y;
    }

    template<multiplicative_monoid M,std::size_t Rank>
    void multiply(tensor_view<M,Rank> &c,const tensor_view<M,Rank> &a,const tensor_view<M,Rank> &b) {
        for (auto && [z,x,y] : zip(c,a,b))
            z = x * y;
    }

    template<multiplicative_group M,std::size_t Rank>
    void divide(tensor_view<M,Rank> &c,const tensor_view<M,Rank> &a,const tensor_view<M,Rank> &b) {
        for (auto && [z,x,y] : zip(c,a,b))
            z = x / y;
    }


    template<additive_monoid M,std::size_t Rank>
    tensor_view<M,Rank> & inplace_add(tensor_view<M,Rank> &b,const tensor_view<M,Rank> &a)
    {
        for (auto && [y,x] : zip(b,a))
            y+=x;
        return b;
    }

    template<additive_group M,std::size_t Rank>
    tensor_view<M,Rank> & inplace_subtract(tensor_view<M,Rank> &b,const tensor_view<M,Rank> &a)
    {
        for (auto && [y,x] : zip(b,a))
            y-=x;
        return b;
    }

    template<multiplicative_monoid M,std::size_t Rank>
    tensor_view<M,Rank> & inplace_multiply(tensor_view<M,Rank> &b,const tensor_view<M,Rank> &a)
    {
        for (auto && [y,x] : zip(b,a))
            y*=x;
        return b;
    }

    template<multiplicative_group M,std::size_t Rank>
    tensor_view<M,Rank> & inplace_divide(tensor_view<M,Rank> &b,const tensor_view<M,Rank> &a)
    {
        for (auto && [y,x] : zip(b,a))
            y/=x;
        return b;
    }

    template<additive_monoid M,std::size_t Rank,std::convertible_to<M> O>
    tensor_view<M,Rank> & inplace_add(tensor_view<M,Rank> &b,const O &k)
    {
        for (auto & y : b)
            y+=k;
        return b;
    }

    template<additive_group M,std::size_t Rank,std::convertible_to<M> O>
    tensor_view<M,Rank> & inplace_subtract(tensor_view<M,Rank> &b,const O &k)
    {
        for (auto & y : b)
            y-=k;
        return b;
    }

    template<multiplicative_monoid M,std::size_t Rank,std::convertible_to<M> O>
    tensor_view<M,Rank> & inplace_multiply(tensor_view<M,Rank> &b,const O &k)
    {
        for (auto &y:b)
            y*=k;
        return b;
    }


    template<multiplicative_group M,std::size_t Rank,std::convertible_to<M> O>
    tensor_view<M,Rank> & inplace_divide(tensor_view<M,Rank> &b,const O &k)
    {
        for (auto &y : b)
            y/=k;
        return b;
    }

    namespace tensor_algebra
    {
        template<additive_monoid M,std::size_t Rank>
        tensor_view<M,Rank>& operator+=(tensor_view<M,Rank>&x,const tensor_view<M,Rank>&y)
        {
            return inplace_add(x,y);
        }

        template<additive_group M,std::size_t Rank>
        tensor_view<M,Rank>& operator-=(tensor_view<M,Rank>&x,const tensor_view<M,Rank>&y)
        {
            return inplace_subtract(x,y);
        }

        template<multiplicative_monoid M,std::size_t Rank>
        tensor_view<M,Rank>& operator*=(tensor_view<M,Rank>&x,const tensor_view<M,Rank>&y)
        {
            return inplace_multiply(x,y);
        }

        template<multiplicative_group M,std::size_t Rank>
        tensor_view<M,Rank>& operator/=(tensor_view<M,Rank>&x,const tensor_view<M,Rank>&y)
        {
            return inplace_multiply(x,y);
        }


        template<additive_monoid M,std::size_t Rank>
        flat_tensor<M,Rank> operator+(const tensor_view<M,Rank>&x,const tensor_view<M,Rank>&y)
        {
            flat_tensor<M,Rank> z(x.shape());
            add(z,x,y);
            return z;
        }

        template<additive_group M,std::size_t Rank>
        flat_tensor<M,Rank> operator-(const tensor_view<M,Rank>&x,const tensor_view<M,Rank>&y)
        {
            flat_tensor<M,Rank> z(x.shape());
            subtract(z,x,y);
            return z;
        }

        template<multiplicative_monoid M,std::size_t Rank>
        flat_tensor<M,Rank> operator*(const tensor_view<M,Rank>&x,const tensor_view<M,Rank>&y)
        {
            flat_tensor<M,Rank> z(x.shape());
            multiply(z,x,y);
            return z;
        }

        template<multiplicative_group M,std::size_t Rank>
        flat_tensor<M,Rank> operator/(const tensor_view<M,Rank>&x,const tensor_view<M,Rank>&y)
        {
            flat_tensor<M,Rank> z(x.shape());
            divide(z,x,y);
            return z;
        }

        template<multiplicative_monoid M,std::size_t Rank, std::convertible_to<M> O>
        flat_tensor<M,Rank> operator*(const tensor_view<M,Rank>&x,const O&k)
        {
            flat_tensor<M,Rank> z(x);
            inplace_multiply(z,k);
            return z;
        }

        template<multiplicative_monoid M,std::size_t Rank, std::convertible_to<M> O>
        flat_tensor<M,Rank> operator*(const O&k,const tensor_view<M,Rank>&x)
        {
            flat_tensor<M,Rank> z(x);
            inplace_multiply(z,k);
            return z;
        }

        template<multiplicative_group M,std::size_t Rank, std::convertible_to<M> O>
        flat_tensor<M,Rank> operator/(const tensor_view<M,Rank>&x, const O&k)
        {
            flat_tensor<M,Rank> z(x);
            inplace_divide(z,k);
            return z;
        }

        template<multiplicative_monoid M,std::size_t Rank, std::convertible_to<M> O>
        tensor_view<M,Rank>& operator*=(tensor_view<M,Rank>&x,const O&k)
        {
            return inplace_multiply(x,k);
        }

        template<multiplicative_group M,std::size_t Rank, std::convertible_to<M> O>
        tensor_view<M,Rank>& operator/=(tensor_view<M,Rank>&x, const O&k)
        {
            return inplace_divide(x,k);
        }
    }
}

#endif //CPLIBRARY_VIEW_OPERATIONS_H

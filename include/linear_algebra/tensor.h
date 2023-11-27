//
// Created by ramizouari on 27/11/23.
//

#ifndef CPLIBRARY_TENSOR_H
#define CPLIBRARY_TENSOR_H
#include <array>
#include <vector>
#include <numeric>
#include "view.h"
namespace cp::linalg
{
    template<typename T,std::size_t n>
    struct tensor_t
    {
        using tensor=std::vector<typename tensor_t<T,n-1>::tensor>;
        std::vector<tensor_t<T,n-1>> U;
        T operator[](const std::array<T,n> &I) const
        {
            std::array<T,n-1> subI={};
            for(int i=1;i<n;i++)
                subI[i-1]=I[i];
            return U[subI];
        }
        explicit operator std::vector<tensor_t<T,n-1>>&() const
        {
            return U;
        }
    };

    template<typename T>
    struct tensor_t<T,0>
    {
        using tensor=T;
        tensor U;
        T operator[](const std::array<T,0>&)
        {
            return U;
        }
        operator const T&() const
        {
            return U;
        }
    };

    template<typename T,std::size_t n>
    using tensor=typename tensor_t<T,n>::tensor;

    template<typename T,std::size_t n>
    T get(const tensor<T,n> &A,std::array<std::size_t,n> I)
    {
        if constexpr (n==0)
            return A;
        else
        {
            std::array<std::size_t, n - 1> subI={};
            for(int i=1;i<n;i++)
                subI[i-1]=I[i];
            return get<T,n-1>(A[I[0]],subI);
        }
    }
    template<typename T,std::size_t n>
    tensor<T,n> reshape(const std::vector<T> &A,std::array<std::size_t,n> shape)
    {
        if constexpr (n==0)
            return A[0];
        else
        {
            int m=A.size()/shape[0];
            std::vector<std::vector<T>> B(shape[0],std::vector<T>(m));
            for(int i=0;i<shape[0];i++)
                for(int j=0;j<m;j++)
                    B[i][j]=A[i*m+j];
            tensor<T,n> R(shape[0]);
            std::array<std::size_t,n-1> subshape={};
            for(int i=1;i<n;i++)
                subshape[i-1]=shape[i];
            for(int i=0;i<shape[0];i++)
                R[i]=reshape<T,n-1>(B[i],subshape);
            return R;
        }
    }

    template<typename R,std::size_t Rank>
    struct flat_tensor : public tensor_view<R,Rank>
    {
        std::vector<R> m_data;
        std::array<std::size_t,Rank> m_shape;
        explicit flat_tensor(std::array<std::size_t,Rank> shape): m_shape(shape)
        {
            m_data.resize(std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<>()));
        }
        R& at(std::array<std::size_t,Rank> indexes) override
        {
            int k=0;
            for(int i=0;i<Rank;i++)
                k= k * m_shape[i] + indexes[i];
            return m_data[k];
        }

        const R& at(std::array<std::size_t,Rank> indexes) const override
        {
            int k=0;
            for(int i=0;i<Rank;i++)
                k= k * m_shape[i] + indexes[i];
            return m_data[k];
        }


        R& at(size_t k)
        {
            return m_data[k];
        }

        const R& at(size_t k) const
        {
            return m_data[k];
        }

        std::size_t size() const override
        {
            return m_data.size();
        }

        const R* data() const
        {
            return m_data.data();
        }

        R* data()
        {
            return m_data.data();
        }

        std::array<std::size_t,Rank> shape() const override
        {
            return m_shape;
        }
    };

    template<typename R>
    struct flat_tensor<R,dynamic_extent> : public tensor_view<R,dynamic_extent>
    {
        std::vector<R> m_data;
        std::vector<std::size_t> m_shape;
        explicit flat_tensor(std::vector<std::size_t> shape):m_shape(shape)
        {
            m_data.resize(std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<>()));
        }
        R& at(std::vector<std::size_t> indexes) override
        {
            std::size_t k=0;
            for(int i=0;i<indexes.size();i++)
                k=k*m_shape[i]+indexes[i];
            return m_data[k];
        }

        const R& at(std::vector<std::size_t> indexes) const override
        {
            std::size_t k=0;
            for(int i=0;i<indexes.size();i++)
                k=k*m_shape[i]+indexes[i];
            return m_data[k];
        }

        R& at(size_t k)
        {
            return m_data[k];
        }

        const R& at(size_t k) const
        {
            return m_data[k];
        }

        std::size_t size() const override
        {
            return m_data.size();
        }

        const R* data() const
        {
            return m_data.data();
        }

        R* data()
        {
            return m_data.data();
        }

        std::vector<std::size_t> shape() const override
        {
            return m_shape;
        }
    };
}

#endif //CPLIBRARY_TENSOR_H

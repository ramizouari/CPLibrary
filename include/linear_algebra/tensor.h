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
    template<typename R,std::size_t Rank>
    struct flat_tensor : tensor_view<R,Rank>
    {
        std::vector<R> m_data;
        std::array<std::size_t,Rank> m_shape;
        explicit flat_tensor(std::array<std::size_t,Rank> shape): m_shape(shape)
        {
            m_data.resize(std::accumulate(shape.begin(), shape.end(), 1, std::multiplies<>()));
        }

        explicit flat_tensor(const tensor_view<R,Rank> & other): m_data(other.begin(),other.end()),m_shape(other.shape())
        {
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

        explicit flat_tensor(const tensor_view<R,dynamic_extent> & other): m_shape(other.shape()),m_data(other.begin(),other.end())
        {
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

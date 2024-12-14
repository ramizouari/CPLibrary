//
// Created by ramizouari on 13/12/24.
//

#ifndef CPLIBRARY_VECTOR_VIEW_H
#define CPLIBRARY_VECTOR_VIEW_H
#include "../view.h"
namespace cp::linalg
{
    template<typename R>
    struct vector_view : tensor_view<R,1>
    {
        R* m_data;
        std::size_t m_size;
        vector_view(R* _data,std::size_t _size):m_data(_data),m_size(_size){}
        vector_view(std::vector<R> &v):m_data(v.data()),m_size(v.size()){}
        template<std::size_t N>
        vector_view(std::array<R,N> &v):m_data(v.data()),m_size(v.size()){}
        ~vector_view() override = default;
        std::size_t size() const override
        {
            return m_size;
        }
        virtual R& at(std::size_t i)
        {
            return m_data[i];
        }

        virtual const R& at(std::size_t i) const
        {
            return m_data[i];
        }
        R& at(std::array<std::size_t,1> indexes) override
        {
            return at(indexes[0]);
        }
        const R& at(std::array<std::size_t,1> indexes) const override
        {
            return at(indexes[0]);
        }
        std::array<std::size_t,1> shape() const override
        {
            return {size()};
        }
        tensor_subview<R,1> slice(std::array<std::size_t,1> start,std::array<std::size_t,1> end) override
        {
            return tensor_subview<R,1>(*this,start,end);
        }

        tensor_subview<R,1> slice(std::array<std::size_t,1> start,std::array<std::size_t,1> end, std::array<std::size_t,1> step) override
        {
            return tensor_subview<R,1>(*this,start,end,step);
        }

        vector_view& operator=(const std::vector<R>& O)
        {
            for(int i=0;i<size();i++)
                at(i)=O.at(i);
            return *this;
        }
    };
}

#endif //CPLIBRARY_VECTOR_VIEW_H

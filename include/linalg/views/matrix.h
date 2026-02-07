//
// Created by ramizouari on 13/12/24.
//

#ifndef CPLIBRARY_MATRIX_VIEW_H
#define CPLIBRARY_MATRIX_VIEW_H


#include <memory>

#include "../../tensors/view.h"
namespace cp::linalg
{
    template<typename R>
    struct matrix_view : tensor_view<R,2>
    {
        std::unique_ptr<R*[]> m_data;
        std::size_t m_rows,m_cols;

        matrix_view(std::size_t rows,std::size_t cols): m_data(std::make_unique<R*[]>(rows)),m_rows(rows),m_cols(cols){}

        matrix_view(R** _data,std::size_t rows,std::size_t cols):matrix_view(rows,cols)
        {
            std::copy(_data,_data+rows,m_data.get());
        }

        matrix_view(std::vector<std::vector<R>> &v):matrix_view(v.size(),v.size()?v.front().size():0)
        {
            for (int i=0;i<m_rows;i++)
                m_data[i] = v[i].data();
        }
        template<std::size_t N>
        matrix_view(std::array<R,N> &v):matrix_view(v.size(),v.size()?v.front().size():0)
        {
            for (int i=0;i<m_rows;i++)
                m_data[i] = v[i].data();
        }
        ~matrix_view() override = default;
        std::size_t size() const override
        {
            return m_rows * m_cols;
        }
        virtual R& at(std::size_t i,std::size_t j)
        {
            return m_data[i][j];
        }
        virtual const R& at(std::size_t i,std::size_t j) const
        {
            return m_data[i][j];
        }
        R& at(std::array<std::size_t,2> indexes) override
        {
            return at(indexes[0],indexes[1]);
        }
        const R& at(std::array<std::size_t,2> indexes) const override
        {
            return at(indexes[0],indexes[1]);
        }
        std::array<std::size_t,2> shape() const override
        {
            return {m_rows,m_cols};
        }
        tensor_subview<R,2> slice(std::array<std::size_t,2> start,std::array<std::size_t,2> end) override
        {
            return tensor_subview<R,2>(*this,start,end);
        }

        tensor_subview<R,2> slice(std::array<std::size_t,2> start,std::array<std::size_t,2> end, std::array<std::size_t,2> step) override
        {
            return tensor_subview<R,2>(*this,start,end,step);
        }

    };
}

#endif //CPLIBRARY_MATRIX_VIEW_H

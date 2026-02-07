//
// Created by ramizouari on 14/12/24.
//

#ifndef FLAT_TENSOR_VIEW_H
#define FLAT_TENSOR_VIEW_H
#include "../../tensors/view.h"

namespace cp::tensors {
    template<typename R,std::size_t Rank>
    struct flat_tensor_view : tensor_view<R,Rank>
    {
        using index_array = tensor_view<R,Rank>::index_array;
        R* m_data;
        index_array m_shape;

        explicit flat_tensor_view(R* m_data,index_array shape): m_data(m_data),m_shape(shape)
        {
        }

        R& at(index_array indexes) override
        {
            int k=0;
            for(int i=0;i<indexes.size();i++) k= k * m_shape[i] + indexes[i];
            return m_data[k];
        }

        const R& at(index_array indexes) const override
        {
            int k=0;
            for(int i=0;i<indexes.size();i++) k= k * m_shape[i] + indexes[i];
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

        const R* data() const
        {
            return m_data;
        }

        R* data()
        {
            return m_data;
        }

        index_array shape() const override
        {
            return m_shape;
        }
    };
}

#endif //FLAT_TENSOR_VIEW_H

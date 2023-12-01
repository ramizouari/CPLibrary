//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_WRAPPERS_H
#define CPLIBRARY_WRAPPERS_H

#include "view.h"

namespace cp::linalg
{
    template<typename T, std::size_t Rank>
    struct zero_padding : public tensor_view<T, Rank>
    {
        const T zero{};
        tensor_view<T, Rank> &src;
        std::array<std::size_t, Rank> sizes;
        zero_padding(tensor_view<T, Rank> &_src, std::array<std::size_t, Rank> _sizes) : src(_src), sizes(_sizes)
        {
        }

        T& at(std::array<std::size_t, Rank> I) override
        {
            for (int i = 0; i < Rank; i++)
                if (I[i] >= sizes[i])
                    return zero;
            return src.at(I);
        }

        const T& at(std::array<std::size_t, Rank> I) const override
        {
            for (int i = 0; i < Rank; i++)
                if (I[i] >= sizes[i])
                    return zero;
            return src.at(I);
        }

        std::array<std::size_t, Rank> shape() const override
        {
            return sizes;
        }

    };
}

#endif //CPLIBRARY_WRAPPERS_H

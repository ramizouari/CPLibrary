//
// Created by ramizouari on 14/12/24.
//

#ifndef CPLIBRARY_LINALG_UTILS_H
#define CPLIBRARY_LINALG_UTILS_H
#include <array>
#include <vector>
#include <ranges>
namespace cp
{

    inline constexpr struct size_tag_t {} size_tag;
    inline constexpr std::size_t dynamic_extent = -1;
    template<typename Range>
    concept sized_random_access = std::ranges::random_access_range<Range> && std::ranges::sized_range<Range>;


    template<sized_random_access P,sized_random_access Q>
    P & increment(P &X, const Q &shape)
    {
        for(int i=0;i<shape.size();i++)
        {
            ++X[i];
            if(X[i]<shape[i])
                break;
            X[i]=0;
        }
        return X;
    }

    template<sized_random_access P,sized_random_access Q>
    bool rincrement(P &X, const Q &shape)
    {
        for(int i=shape.size()-1;i>=0;--i)
        {
            ++X[i];
            if(X[i]<shape[i]) return false;
            X[i]=0;
        }
        return true;
    }

    template<sized_random_access P,sized_random_access Q>
    P & decrement(P &X, const Q &shape)
    {
        for(int i=0;i<shape.size();--i)
        {
            if (X[i] > 0) {
                --X[i];
                break;
            }
            X[i]=shape[i]-1;
        }
        return X;
    }


    template<sized_random_access P,sized_random_access Q>
    P & rdecrement(P &X, const Q &shape)
    {
        for(int i=shape.size()-1;i>=0;--i)
        {
            if (X[i] > 0) {
                --X[i];
                break;
            }
            X[i]=shape[i]-1;
        }
        return X;
    }

}
#endif //CPLIBRARY_LINALG_UTILS_H

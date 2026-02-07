//
// Created by ramizouari on 14/12/24.
//

#ifndef CPLIBRARY_LINALG_UTILS_H
#define CPLIBRARY_LINALG_UTILS_H
#include <generator>


namespace cp
{

    inline constexpr struct size_tag_t {} size_tag;
    inline constexpr std::size_t dynamic_extent = -1;


    namespace linalg {
        using cp::dynamic_extent;
    }

    template<std::size_t... ext>
    concept all_dynamic = ((ext == dynamic_extent) && ... );

    template<std::size_t... ext>
    concept none_dynamic = ((ext != dynamic_extent) && ... );

}
#endif //CPLIBRARY_LINALG_UTILS_H

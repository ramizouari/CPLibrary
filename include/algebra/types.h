//
// Created by ramizouari on 13/12/24.
//

#ifndef ALGEBRAIC_TYPES_H
#define ALGEBRAIC_TYPES_H
#include <cstdint>

namespace cp
{
    using natural = std::uint64_t;
    using integer = std::int64_t;
    using real = long double;
    using IR=real;
    using IC= std::complex<IR>;
}

#endif //ALGEBRAIC_TYPES_H

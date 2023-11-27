//
// Created by ramizouari on 22/11/23.
//

#ifndef CPLIBRARY_VARINT_H
#define CPLIBRARY_VARINT_H
#include <cstdint>
#include <array>
#include <vector>
#include <istream>
#include <ostream>
#include <algorithm>
#include <bitset>
#include <iomanip>
#include "utils.h"

namespace cp
{
    struct varint
    {
        void reduce();
        std::vector<std::uint64_t> A;
        varint() = default;
        varint(std::uint64_t x);
        varint operator+(const varint &B) const;
        varint& operator+=(const varint &B);
        varint operator-(const varint &B) const;
        varint& operator-=(const varint &B);
        varint operator*(const varint &B) const;
        varint& operator*=(const varint &B);
        varint operator<<(int k) const;
        varint operator>>(int k) const;
        varint& operator<<=(unsigned k);
        varint& operator>>=(unsigned k);
        [[nodiscard]] std::pair<varint,varint> euclidean_division(const varint& B) const;
        varint operator/(const varint &B) const;
        varint operator%(const varint &B) const;
        varint& operator/=(const varint &B);
        varint& operator%=(const varint &B);
        std::strong_ordering operator<=>(const varint &B) const;
        varint operator~() const;
        varint operator-() const;
        varint& operator++();
        varint operator++(int);
        varint& operator--();
        varint operator--(int);
        varint& operator&=(const varint &B);
        varint& operator|=(const varint &B);
        varint& operator^=(const varint &B);
        varint operator&(const varint &B) const;
        varint operator|(const varint &B) const;
        varint operator^(const varint &B) const;
        explicit operator std::uint64_t () const;
        bool operator==(const varint &B) const= default;
    };
    std::ostream& operator<<(std::ostream &H, varint A);
    std::istream& operator>>(std::istream &H, varint &A);
}
#endif //CPLIBRARY_VARINT_H

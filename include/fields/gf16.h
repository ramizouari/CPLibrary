//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_GF16_H
#define CPLIBRARY_GF16_H

#include <cstdint>

#include "algebra/abstract_algebra.h"
#include "utils.h"
namespace cp
{
    struct GF16
    {
        std::uint8_t v;
        GF16(integer a=0,integer b=0,integer c=0,integer d=0):v((a&1)|((b&1)<<1)|((c&1)<<2)|((d&1)<<3)){}
        GF16 operator+(const GF16 &a) const
        {
            return GF16(representation(v^a.v));
        }
        GF16 operator*(const GF16 &a) const
        {
            auto z=xor_multiply(v,a.v);
            return GF16(representation((z>>4) | (z>>3) | z));
        }
        GF16 operator-() const
        {
            return GF16(representation(v));
        }
        GF16 inv() const
        {
            return GF16(representation(v^1));
        }
        bool operator==(const GF16 &a) const = default;
        GF16& operator+=(const GF16 &a)
        {
            v^=a.v;
            return *this;
        }
        GF16& operator*=(const GF16 &a)
        {
            v=(v&a.v)^((v^a.v)&1);
            return *this;
        }
        GF16& operator-=(const GF16 &a)
        {
            v^=a.v;
            return *this;
        }
        GF16& operator/=(const GF16 &a)
        {

            return *this;
        }

        GF16 operator/(const GF16 &a) const
        {
            return GF16(representation(v^a.v));
        }

    protected:
        struct representation
        {
            std::uint8_t v;
        };
        GF16(representation r):v(r.v){}
    };
}

#endif //CPLIBRARY_GF16_H

//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_GF8_H
#define CPLIBRARY_GF8_H

#include "algebra/abstract_algebra.h"
#include "utils.h"

namespace cp
{
    struct GF8
    {
        std::uint8_t v;
        GF8(integer a=0,integer b=0,integer c=0):v((a&1)|((b&1)<<1)|((c&1)<<2)){}
        GF8 operator+(const GF8 &a) const
        {
            return GF8(representation(v^a.v));
        }
        GF8 operator*(const GF8 &a) const
        {
            auto z=xor_multiply(v,a.v);
            return GF8(representation((z>>4) | (z>>3) | z));
        }
        GF8 operator-() const
        {
            return GF8(representation(v));
        }
        GF8 inv() const
        {
            return GF8(representation(v^1));
        }
        bool operator==(const GF8 &a) const = default;
        GF8& operator+=(const GF8 &a)
        {
            v^=a.v;
            return *this;
        }
        GF8& operator*=(const GF8 &a)
        {
            v=(v&a.v)^((v^a.v)&1);
            return *this;
        }
        GF8& operator-=(const GF8 &a)
        {
            v^=a.v;
            return *this;
        }
        GF8& operator/=(const GF8 &a)
        {

            return *this;
        }

        GF8 operator/(const GF8 &a) const
        {
            return GF8(representation(v^a.v));
        }
    protected:
        struct representation
        {
            std::uint8_t z;
        };
        GF8(representation r):v(r.z){}
    };
}

#endif //CPLIBRARY_GF8_H

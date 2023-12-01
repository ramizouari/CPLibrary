//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_GF4_H
#define CPLIBRARY_GF4_H

#include "algebra/abstract_algebra.h"
#include "bit"
namespace cp
{
    struct GF4
    {
        std::uint8_t v;
        GF4(integer a=0,integer b=0):v((a&1)|((b&1)<<1)){}
        GF4 operator+(const GF4 &a) const
        {
            return GF4(representation(v^a.v));
        }
        GF4 operator*(const GF4 &a) const
        {
            auto p= v && a.v ? ((a.v+v - 2) %3)+1 :0;
            return GF4(representation(p));
        }
        GF4 operator-() const
        {
            return GF4(representation(v));
        }
        GF4 inv() const
        {
            return GF4(representation(v^1));
        }
        bool operator==(const GF4 &a) const = default;
        GF4& operator+=(const GF4 &a)
        {
            v^=a.v;
            return *this;
        }
        GF4& operator*=(const GF4 &a)
        {
            v=((*this)*a).v;
            return *this;
        }
        GF4& operator-=(const GF4 &a)
        {
            v^=a.v;
            return *this;
        }
        GF4& operator/=(const GF4 &a)
        {

            return *this;
        }

        GF4 operator/(const GF4 &a) const
        {
            auto p=(v&1)?a.v:0;
            auto q=(v&2) && (a.v&2)?0b11:0;
            auto r=(v&2) && (a.v&1)?0b10:0;
            return GF4(representation(p^q^r));        }
    protected:
        struct representation
        {
            std::uint8_t v;
        };
        GF4(representation r):v(r.v){}
    };
}

#endif //CPLIBRARY_GF4_H

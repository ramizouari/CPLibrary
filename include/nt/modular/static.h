//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_STATIC_H
#define CPLIBRARY_STATIC_H
#include "fixed.h"
namespace cp
{
    inline static constexpr integer static_modulus=-1;
    template<>
    struct cyclic<static_modulus>
    {
        integer n;
        inline static integer m=1;
        inline static bool assume_prime=true;
        constexpr cyclic(integer o=0):n((o+m)%m){}
        bool operator==(integer O) const
        {
            return n==(m+O)%m;
        }

        bool operator==(const cyclic &O) const = default;

        cyclic& operator+=(const cyclic &O)
        {
            n=(n+O.n)%m;
            return *this;
        }
        cyclic& operator-=(const cyclic &O)
        {
            n=(n+m-O.n)%m;
            return *this;
        }

        cyclic& operator*=(const cyclic &O)
        {
            n=(n*O.n)%m;
            return *this;
        }

        cyclic inv() const
        {
            if(assume_prime)
                return pow(*this,m-2);
            else return pinv();
        }

        cyclic& operator/=(const cyclic &O)
        {
            return (*this)*=O.inv();
        }

        cyclic operator-() const
        {
            return m-n;
        }

        cyclic pinv() const
        {
            return egcd(n,m).a;
        }

        explicit operator integer&()
        {
            return n;
        }

        explicit operator const integer&() const
        {
            return n;
        }

        static integer modulus()
        {
            return m;
        }
        static void set_modulus(integer _m)
        {
            m=_m;
        }
    };

    using d_cyclic=cyclic<static_modulus>;

}

#endif //CPLIBRARY_STATIC_H

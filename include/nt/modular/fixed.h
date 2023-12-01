//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_FIXED_H
#define CPLIBRARY_FIXED_H
#include "algebra/abstract_algebra.h"


namespace cp
{
    template<integer mod>
    struct cyclic
    {
        integer n;
        inline static bool assume_prime=true;
        inline static constexpr integer m = mod;
        constexpr cyclic(integer o=0):n((o+m)%m){}
        bool operator==(int O) const
        {
            return n==(m+O)%m;
        }

        bool operator!=(int O) const
        {
            return n!=(m+O)%m;
        }

        bool operator==(cyclic O) const
        {
            return n==O.n;
        }

        bool operator!=(cyclic O) const
        {
            return n!=O.n;
        }

        cyclic operator-() const
        {
            return cyclic(-n);
        }

        auto& operator+=(const cyclic &O)
        {
            n=(n+O.n)%m;
            return *this;
        }
        auto& operator-=(const cyclic &O)
        {
            n=(n+m-O.n)%m;
            return *this;
        }

        auto& operator*=(const cyclic &O)
        {
            n=(n*O.n)%m;
            return *this;
        }

        auto& operator/=(const cyclic &O)
        {
            return (*this)*=O.inv();
        }

        auto operator*(const cyclic &O) const
        {
            auto w=*this;
            return w*=O;
        }

        auto operator+(const cyclic &O) const
        {
            auto w=*this;
            return w+=O;
        }

        auto operator-(const cyclic &O) const
        {
            auto w=*this;
            return w-=O;
        }

        auto operator/(const cyclic &O) const
        {
            return (*this)*O.inv();
        }

        cyclic pinv() const
        {
            return egcd(n,m).a;
        }

        auto inv() const
        {
            if(assume_prime)
                return pow(*this,m-2);
            else return pinv();
        }

        auto& operator++()
        {
            return *this+=1;
        }

        auto& operator--()
        {
            return *this-=1;
        }

        auto operator++(int)
        {
            cyclic r(n);
            *this += 1;
            return r;
        }

        auto operator--(int)
        {
            cyclic r(n);
            *this -= 1;
            return r;
        }

        explicit operator integer&()
        {
            return n;
        }

        explicit operator const integer&() const
        {
            return n;
        }

        static constexpr integer modulus()
        {
            return m;
        }
    };

    template<integer m>
    auto operator*(integer k,cyclic<m> s)
    {
        return s*k;
    }

    template<integer m>
    auto operator+(integer k,cyclic<m> s)
    {
        return s+k;
    }

    template<integer m>
    auto operator-(integer k,cyclic<m> s)
    {
        return (-s)+k;
    }
}

#endif //CPLIBRARY_FIXED_H

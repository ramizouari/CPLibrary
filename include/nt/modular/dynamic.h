//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_DYNAMIC_H
#define CPLIBRARY_DYNAMIC_H
#include "fixed.h"
namespace cp
{
    inline static constexpr integer dynamic_modulus=-2;
    template<>
    struct cyclic<dynamic_modulus>
    {
        integer m,n;
        bool assume_prime=true;
        cyclic(integer o=0,integer q=0):m(q),n(m?(o+m)%m:o){}
        bool operator==(integer O) const
        {
            if(!m) return n==O;
            return n==(m+O%m)%m;
        }

        bool operator==(const cyclic &O) const = default;

        cyclic& operator+=(const cyclic &O)
        {
            if(!m) m=O.m;
            n+=O.n;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator+=(integer O)
        {
            n=n+O;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator-=(const cyclic &O)
        {
            if(!m)
                m=O.m;
            n+=m-O.n;
            if(m)
                n%=m;
            return *this;
        }

        cyclic& operator-=(integer O)
        {
            n+=m-O;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator*=(const cyclic &O)
        {
            if(!m) m=O.m;
            n*=O.n;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator*=(integer O)
        {
            n*=O;
            if(m) n%=m;
            return *this;
        }

        cyclic& operator=(integer O)
        {
            n=O;
            if(m) n%=m;
            return *this;
        }

        cyclic inv() const
        {
            if(m==1)
                return *this;
            else if(assume_prime)
                return pow(*this,m-2,m);
            else return pinv();
        }

        cyclic& operator/=(const cyclic &O)
        {
            return (*this)*=O.inv();
        }

        cyclic& operator/=(integer O)
        {
            return (*this)*=cyclic(O,m).inv();
        }

        cyclic operator-() const
        {
            return {m-n,m};
        }

        cyclic pinv() const
        {
            return {egcd(n,m).a,m};
        }

        explicit operator integer&()
        {
            return n;
        }

        explicit operator const integer&() const
        {
            return n;
        }

        integer modulus() const
        {
            return m;
        }
    };

}

#endif //CPLIBRARY_DYNAMIC_H

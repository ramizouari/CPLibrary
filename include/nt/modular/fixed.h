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
        static constexpr integer m = mod;
        constexpr cyclic(integer o=0):n((o%m+m)%m){}
        bool operator==(integer O) const
        {
            return n==(m+O%m)%m;
        }

        bool operator==(const cyclic &O) const = default;

        cyclic operator-() const
        {
            return cyclic(m-n);
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

        cyclic pinv() const
        {
            return egcd(n,m).a;
        }

        auto inv() const
        {
            if(assume_prime) return pow(*this,m-2);
            return pinv();
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
    using cyc=cyclic<m>;

    template<typename Cyc,integer m>
    concept ToCyclic=std::convertible_to<Cyc,cyc<m>>;

    template<typename Cyc,integer m>
    concept ToCyclicProper=std::convertible_to<Cyc,cyc<m>> && !std::same_as<Cyc,cyc<m>>;

    template<integer m,ToCyclic<m> O>
    cyc<m> operator*(const cyc<m> &A,const O &B)
    {
        auto C=A;
        return C*=B;
    }

    template<integer m,ToCyclicProper<m> O>
    cyc<m> operator*(const O & A,const cyc<m> & B)
    {
        cyc<m> C=A;
        return C*=B;
    }

    template<integer m,ToCyclic<m> O>
    cyc<m> operator+(const cyc<m> &A,const O &B)
    {
        cyc<m> C=A;
        return C+=B;
    }

    template<integer m,ToCyclicProper<m> O>
    cyc<m> operator+(const O &A,const cyc<m> &B)
    {
        cyc<m> C=A;
        return C+=B;
    }

    template<integer m,ToCyclic<m> O>
    cyc<m> operator-(const cyc<m>& A,const O &B)
    {
        cyc<m> C=A;
        return C-=B;
    }

    template<integer m,ToCyclicProper<m> O>
    cyc<m> operator-(const O& A,const cyc<m> &B)
    {
        cyc<m> C=A;
        return C-=B;
    }


    template<integer m,ToCyclic<m> O>
    cyc<m> operator/(const cyclic<m> &A,const O &B)
    {
        cyc<m> C=A;
        return C/=B;
    }

    template<integer m,ToCyclicProper<m> O>
    cyc<m> operator/(const O &A,const cyc<m> &B)
    {
        cyc<m> C=A;
        return C/=B;
    }

    template<integer m>
    cyc<m>& operator++(cyc<m> & x)
    {
        return x+=1;
    }


    template<integer m>
    cyc<m>& operator--(cyc<m> & x)
    {
        return x-=1;
    }

    template<integer m>
    cyc<m> operator++(cyc<m>&x,int)
    {
        auto y=x;
        x += 1;
        return y;
    }

    template<integer m>
    cyc<m> operator--(cyc<m>&x,int)
    {
        auto y=x;
        x -= 1;
        return y;
    }

    template<integer m>
    std::ostream& operator<<(std::ostream & H,const cyc<m> &x) {
        return H << x.n;
    }

    template<integer m>
    std::istream& operator>>(std::istream & H, cyc<m> &x) {
        integer n;
        H >> n;
        x=n;
        return H;
    }
}

#endif //CPLIBRARY_FIXED_H

//
// Created by ramizouari on 22/11/23.
//

#ifndef CPLIBRARY_OPERATIONS_H
#define CPLIBRARY_OPERATIONS_H
#include "varint.h"
#include "bigint.h"

namespace cp
{
    template<typename T,std::size_t N,std::size_t R>
    using remove_cvref_if_less = std::conditional_t<N<R,std::remove_cvref_t<T>,T>;

    template<std::size_t N,std::size_t M>
    std::strong_ordering operator<=>(const big_int<N> &A,const big_int<M> &B)
    {
        using T1=remove_cvref_if_less<const big_int<N>&,N,M>;
        using T2=remove_cvref_if_less<const big_int<M>&,M,N>;
        return static_cast<T1>(A)<=>static_cast<T2>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<std::max(N,M)> operator+(const big_int<N> &A,const big_int<M> &B)
    {
        if constexpr (N < M)
            return B+A;
        else
            return A+static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<N>& operator+=(big_int<N> &A,const big_int<M> &B)
    {
        return A+=static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<std::max(N,M)> operator-(const big_int<N> &A,const big_int<M> &B)
    {
        if constexpr (N < M)
            return -B+A;
        else
            return A-static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<N>& operator-=(big_int<N> &A,const big_int<M> &B)
    {
        return A-=static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<std::max(N,M)> operator*(const big_int<N> &A,const big_int<M> &B)
    {
        if constexpr (N < M)
            return B*A;
        else
            return A*static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<N>& operator*=(big_int<N> &A,const big_int<M> &B)
    {
        return A*=static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<std::max(N,M)> operator/(const big_int<N> &A,const big_int<M> &B)
    {
        constexpr std::size_t R=std::max(N,M);
        using T1=remove_cvref_if_less<const big_int<R>&,N,R>;
        using T2=remove_cvref_if_less<const big_int<R>&,M,R>;
        return static_cast<T1>(A)/static_cast<T2>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<N>& operator/=(big_int<N> &A,const big_int<M> &B)
    {
        return A/=static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<std::max(N,M)> operator%(const big_int<N> &A,const big_int<M> &B)
    {
        constexpr std::size_t R=std::max(N,M);
        using T1=remove_cvref_if_less<const big_int<R>&,N,R>;
        using T2=remove_cvref_if_less<const big_int<R>&,M,R>;
        return static_cast<T1>(A)%static_cast<T2>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<N>& operator%=(big_int<N> &A,const big_int<M> &B)
    {
        return A%=static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<std::max(N,M)> operator&(const big_int<N> &A,const big_int<M> &B)
    {
        if constexpr (N < M)
            return B&A;
        else
            return A&static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<std::max(N,M)> operator|(const big_int<N> &A,const big_int<M> &B)
    {
        if constexpr (N < M)
            return B|A;
        else
            return A|static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<N>& operator&=(big_int<N> &A,const big_int<M> &B)
    {
        return A&=static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<N>& operator|=(big_int<N> &A,const big_int<M> &B)
    {
        return A|=static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<std::max(N,M)> operator^(const big_int<N> &A,const big_int<M> &B)
    {
        if constexpr (N < M)
            return B^A;
        else
            return A^static_cast<big_int<N>>(B);
    }

    template<std::size_t N,std::size_t M>
    big_int<N>& operator^=(big_int<N> &A,const big_int<M> &B)
    {
        return A^=static_cast<big_int<N>>(B);
    }

    template<std::size_t N>
    varint operator+(const varint &A,const big_int<N> &B)
    {
        varint C;
        std::copy(B.A.begin(),B.A.end(),std::back_inserter(C.A));
        return A+C;
    }

    template<std::size_t N>
    varint operator+(const big_int<N> &A,const varint &B)
    {
        return B+A;
    }

    template<std::size_t N>
    varint operator-(const varint &A,const big_int<N> &B)
    {
        varint C;
        std::copy(B.A.begin(),B.A.end(),std::back_inserter(C.A));
        return A-C;
    }

    template<std::size_t N>
    varint operator-(const big_int<N> &A,const varint &B)
    {
        varint C;
        std::copy(A.A.begin(),A.A.end(),std::back_inserter(C.A));
        return C-B;
    }

    template<std::size_t N>
    varint operator*(const varint &A,const big_int<N> &B)
    {
        varint C;
        std::copy(B.A.begin(),B.A.end(),std::back_inserter(C.A));
        return A*C;
    }

    template<std::size_t N>
    varint operator*(const big_int<N> &A,const varint &B)
    {
        return B*A;
    }

    template<std::size_t N>
    varint operator/(const varint &A,const big_int<N> &B)
    {
        varint C;
        std::copy(B.A.begin(),B.A.end(),std::back_inserter(C.A));
        return A/C;
    }

    template<std::size_t N>
    varint operator/(const big_int<N> &A,const varint &B)
    {
        varint C;
        std::copy(A.A.begin(),A.A.end(),std::back_inserter(C.A));
        return C/B;
    }

    template<std::size_t N>
    varint operator%(const varint &A,const big_int<N> &B)
    {
        varint C;
        std::copy(B.A.begin(),B.A.end(),std::back_inserter(C.A));
        return A%C;
    }

    template<std::size_t N>
    varint operator%(const big_int<N> &A,const varint &B)
    {
        varint C;
        std::copy(A.A.begin(),A.A.end(),std::back_inserter(C.A));
        return C%B;
    }

    template<std::size_t N>
    varint operator&(const varint &A,const big_int<N> &B)
    {
        varint C;
        std::copy(B.A.begin(),B.A.end(),std::back_inserter(C.A));
        return A&C;
    }

    template<std::size_t N>
    varint operator&(const big_int<N> &A,const varint &B)
    {
        return B&A;
    }

    template<std::size_t N>
    varint operator|(const varint &A,const big_int<N> &B)
    {
        varint C;
        std::copy(B.A.begin(),B.A.end(),std::back_inserter(C.A));
        return A|C;
    }

    template<std::size_t N>
    varint operator|(const big_int<N> &A,const varint &B)
    {
        return B|A;
    }

    template<std::size_t N>
    varint operator^(const varint &A,const big_int<N> &B)
    {
        varint C;
        std::copy(B.A.begin(),B.A.end(),std::back_inserter(C.A));
        return A^C;
    }

    template<std::size_t N>
    varint operator^(const big_int<N> &A,const varint &B)
    {
        return B^A;
    }

    template<std::size_t N>
    varint& operator+=(varint &A,const big_int<N> &B)
    {
        return A+=static_cast<varint>(B);
    }

    template<std::size_t N>
    varint& operator-=(varint &A,const big_int<N> &B)
    {
        return A-=static_cast<varint>(B);
    }

    template<std::size_t N>
    varint& operator*=(varint &A,const big_int<N> &B)
    {
        return A*=static_cast<varint>(B);
    }

    template<std::size_t N>
    varint& operator/=(varint &A,const big_int<N> &B)
    {
        return A/=static_cast<varint>(B);
    }

    template<std::size_t N>
    varint& operator%=(varint &A,const big_int<N> &B)
    {
        return A%=static_cast<varint>(B);
    }

    template<std::size_t N>
    varint& operator&=(varint &A,const big_int<N> &B)
    {
        return A&=static_cast<varint>(B);
    }

    template<std::size_t N>
    varint& operator|=(varint &A,const big_int<N> &B)
    {
        return A|=static_cast<varint>(B);
    }

    template<std::size_t N>
    varint& operator^=(varint &A,const big_int<N> &B)
    {
        return A^=static_cast<varint>(B);
    }
}


#endif //CPLIBRARY_OPERATIONS_H

//
// Created by ramizouari on 12/12/23.
//

#ifndef CPLIBRARY_BIG_INT_FAST_MULTIPLICATION_H
#define CPLIBRARY_BIG_INT_FAST_MULTIPLICATION_H
#include "bigint.h"
#include "varint.h"
#include "polynomial/fast_polynomial.h"

namespace cp
{
    template<std::size_t N1,std::size_t N2>
    big_int<std::max(N1,N2)> fast_multiplication(const big_int<N1> &A,const big_int<N2> &B)
    {
        constexpr std::size_t N=std::max(N1,N2);
        using H=big_int<N>;
        constexpr integer M=998244353;
        constexpr integer R=H::N;
        std::vector<cyclic<M>> a(4*R),b(4*R);
        for(int i=0;i<4*R;i++)
        {
            a[i]=A.A[i>>2]/(1ull<<(16*(i&3)));
            b[i]=B.A[i>>2]/(1ull<<(16*(i&3)));
        }
        auto c=fast_multiplication(a,b);
        H C;
        for(int i=0;i<4*R;i++)
            C.A[i>>2]+=static_cast<std::uint64_t>(static_cast<integer>(c[i]))<<(16*(i&3));
        return C;
    }

    varint fast_multiplication(const varint &A,const varint &B)
    {
        using H=varint;
        constexpr integer M1=998244353,M2=469762049;
        integer R=std::max(A.A.size(),B.A.size());
        std::vector<cyclic<M1>> a1(4*R),b1(4*R);
        std::vector<cyclic<M2>> a2(4*R),b2(4*R);
        for(int i=0;i<4*R;i++)
        {
            auto r=(16*(i&3));
            auto x=(A.A[i>>2]>>r)&0xFFFF;
            auto y=(B.A[i>>2]>>r)&0xFFFF;
            a1[i]=x;
            b1[i]=y;
            a2[i]=x;
            b2[i]=y;
        }
        auto c1=fast_multiplication(a1,b1);
        auto c2=fast_multiplication(a2,b2);
        std::vector<std::int64_t> c(4*R);
        for(int i=0;i<4*R;i++) {
            auto x=static_cast<integer>(c1[i]);
            auto y=static_cast<integer>(c2[i]);
            //TODO: standardize chinese remainder theorem
            c[i] = chinese_remainder<__int128>({x,y}, {M1, M2});
        }
        H C;
        C.A.resize(4*R+1);
        for(int i=0;i<4*R;i++) for (int j :{0,1,2,3})
        {
            std::uint64_t X=static_cast<std::uint64_t>(static_cast<integer>(c[i])) >> (16*((i+j)&3));
            X&=0xFFFF;
            std::uint64_t carry=0;
            C.A[(i+j)>>2]=add_with_carry(C.A[(i+j)>>2],X << (16*((i+j)&3)),carry);
            if(carry)
                C.A[((i+j)>>2)+1]++;
        }
        C.reduce();
        return C;
    }
}

#endif //CPLIBRARY_BIG_INT_FAST_MULTIPLICATION_H

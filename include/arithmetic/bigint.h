//
// Created by ramizouari on 22/11/23.
//

#ifndef CPLIBRARY_BIGINT_H
#define CPLIBRARY_BIGINT_H
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
    template<std::size_t M>
    struct big_int
    {
        static constexpr std::size_t N=(M+63)/64;
        std::array<std::uint64_t,N> A;
        big_int() = default;
        big_int(std::uint64_t x)
        {
            for(int i=0;i<N;i++)
                A[i]=0;
            A[0]=x;
        }

        template<std::size_t O>
        explicit big_int(const big_int<O> &B)
        {
            for(int i=0;i<N;i++)
                A[i]=0;
            for(int i=0;i<std::min(N,big_int<O>::N);i++)
                A[i]=B.A[i];
        }

        big_int operator+(const big_int &B) const
        {
            big_int C;
            std::uint64_t carry=0;
            for(int i=0;i<N;i++)
            {
                C.A[i]=A[i]+B.A[i]+carry;
                carry=(C.A[i]<A[i])||(C.A[i]<B.A[i]);
            }
            return C;
        }

        big_int& operator+=(const big_int &B)
        {
            std::uint64_t carry=0;
            for(int i=0;i<N;i++)
            {
                A[i]+=B.A[i]+carry;
                carry=(A[i]<B.A[i])||(A[i]<carry);
            }
            return *this;
        }

        big_int& operator-=(const big_int &B)
        {
            std::uint64_t carry=0;
            for(int i=0;i<N;i++)
            {
                A[i]-=B.A[i]+carry;
                carry=(A[i]>B.A[i])||(A[i]>carry);
            }
            return *this;
        }

        big_int operator~() const
        {
            big_int C;
            for(int i=0;i<N;i++)
                C.A[i]=~A[i];
            return C;
        }

        big_int operator-() const
        {
            auto C=~(*this);
            return C+=1;
        }

        big_int operator&(const big_int &B) const
        {
            big_int C;
            for(int i=0;i<N;i++)
                C.A[i]=A[i]&B.A[i];
            return C;
        }

        big_int operator|(const big_int &B) const
        {
            big_int C;
            for(int i=0;i<N;i++)
                C.A[i]=A[i]|B.A[i];
            return C;
        }

        big_int operator^(const big_int &B) const
        {
            big_int C;
            for(int i=0;i<N;i++)
                C.A[i]=A[i]^B.A[i];
            return C;
        }

        big_int& operator&=(const big_int &B)
        {
            for(int i=0;i<N;i++)
                A[i]&=B.A[i];
            return *this;
        }

        big_int& operator|=(const big_int &B)
        {
            for(int i=0;i<N;i++)
                A[i]|=B.A[i];
            return *this;
        }

        big_int& operator^=(const big_int &B)
        {
            for(int i=0;i<N;i++)
                A[i]^=B.A[i];
            return *this;
        }

        big_int operator++(int)
        {
            big_int C=*this;
            *this+=1;
            return C;
        }

        big_int operator--(int)
        {
            big_int C=*this;
            *this-=1;
            return C;
        }

        big_int& operator++()
        {
            return *this+=1;
        }

        big_int& operator--()
        {
            return *this-=1;
        }

        big_int operator-(const big_int &B) const
        {
            return *this+(-B);
        }

        big_int operator*(const big_int &B) const
        {
            big_int C{};
            for(int i=0;i<N;i++)
            {
                std::uint64_t D1=0,D2=0;
                for(int j=0;j<N-i;j++)
                {
                    C.A[i+j]= reduce_carry(C.A[i+j],D1);
                    if(j<N-i-1) C.A[i+j+1]= reduce_carry(C.A[i+j+1],D2);
                    auto Y1=A[i]>>32;
                    auto X1=A[i]&0xffffffff;
                    auto Y2=B.A[j]>>32;
                    auto X2=B.A[j]&0xffffffff;
                    auto Z1=X1*Y2;
                    auto Z2=X2*Y1;
                    C.A[i+j]= add_with_carry(C.A[i+j],X1*X2,D1);
                    C.A[i+j]= add_with_carry(C.A[i+j],Z1 << 32,D1);
                    C.A[i+j]= add_with_carry(C.A[i+j],Z2 << 32,D1);
                    if(j<N-i-1)
                    {
                        C.A[i+j+1] = add_with_carry(C.A[i + j + 1], Y1*Y2, D2);
                        C.A[i + j + 1] = add_with_carry(C.A[i + j + 1], Z1 >> 32 , D2);
                        C.A[i + j + 1] = add_with_carry(C.A[i + j + 1], Z2 >> 32, D2);
                    }
                }
            }
            return C;
        }

        big_int& operator*=(const big_int &B)
        {
            auto X=(*this)*B;
            A=std::move(X.A);
            return *this;
        }

        big_int shift_8r(unsigned r) const
        {
            big_int C;
            for(int i=r;i<N;i++)
                C.A[i-r]=A[i];
            std::fill(C.A.begin()+N-r,C.A.end(),0);
            return C;
        }

        big_int shift_8l(unsigned r) const
        {
            big_int C;
            for(int i=0;i<N-r;i++)
                C.A[i+r]=A[i];
            std::fill(C.A.begin(),C.A.begin()+r,0);
            return C;
        }

        std::pair<big_int,big_int> euclidean_division(const big_int& B) const
        {
            using namespace bits;
            std::bitset<64*N> X,Y;
            for(int i=0;i<N;i++) for(int j=0;j<64;j++)
                {
                    X[i*64+j]=(A[i]>>j)&1;
                    Y[i*64+j]=(B.A[i]>>j)&1;
                }
            std::bitset<64*N> Q,R;
            for(int i=64*N-1;i>=0;i--)
            {
                R<<=1;
                R[0]=X[i];
                if(R>=Y)
                {
                    R-=Y;
                    Q[i]=1;
                }
            }
            big_int C{},D{};
            for(int i=0;i<N;i++) for(int j=0;j<64;j++)
                {
                    C.A[i]|=static_cast<std::uint64_t>(Q[i*64+j])<<j;
                    D.A[i]|=static_cast<std::uint64_t>(R[i*64+j])<<j;
                }
            return {C,D};
        }

        big_int operator/(const big_int &B) const
        {
            auto [X,Y]=euclidean_division(B);
            return X;
        }

        big_int& operator/=(const big_int &B)
        {
            auto X=(*this)/B;
            A=std::move(X.A);
            return *this;
        }

        big_int operator%(const big_int &B) const
        {
            auto [X,Y]=euclidean_division(B);
            return Y;
        }
        big_int & operator%=(const big_int &B)
        {
            auto X=(*this)%B;
            A=std::move(X.A);
            return *this;
        }

        big_int operator<<(int k) const
        {
            big_int C=*this;
            C<<=k;
            return C;
        }

        big_int operator>>(int k) const
        {
            big_int C;
            for(int i=0;i<N;i++)
                C.A[i]=A[i];
            C>>=k;
            return C;
        }

        big_int& operator<<=(int k)
        {
            for(int i=N-1;i>=k;i--)
                A[i]=A[i-k];
            for(int i=0;i<k;i++)
                A[i]=0;
            return *this;
        }

        big_int& operator>>=(int k)
        {
            for(int i=0;i<N-k;i++)
                A[i]=A[i+k];
            for(int i=N-k;i<N;i++)
                A[i]=0;
            return *this;
        }
        std::strong_ordering operator<=>(const big_int &B) const
        {
            return std::lexicographical_compare_three_way(A.rbegin(),A.rend(),B.A.rbegin(),B.A.rend());
        }

        bool operator==(const big_int &B) const = default;

        explicit operator std::uint64_t () const
        {
            return A[0];
        }

        template<std::size_t O>
        explicit operator big_int<O>() const
        {
            return big_int<O>(*this);
        }
    };

    template<std::size_t M>
    big_int<M> operator*(std::uint64_t k,const big_int<M>& A)
    {
        return A*k;
    }

    template<std::size_t M>
    big_int<M> operator+(std::uint64_t k,const big_int<M>& A)
    {
        return A+k;
    }

    template<std::size_t M>
    big_int<M> operator-(std::uint64_t k,const big_int<M>& A)
    {
        return big_int<M>(k)-A;
    }

    template<std::size_t M>
    big_int<M> operator/(std::uint64_t k,const big_int<M>& A)
    {
        return big_int<M>(k)/A;
    }

    template<std::size_t M>
    big_int<M> operator%(std::uint64_t k,const big_int<M>& A)
    {
        return big_int<M>(k)%A;
    }

    template<std::size_t M>
    std::ostream& operator<<(std::ostream &H, big_int<M> A)
    {
        std::uint64_t x=0;
        constexpr std::uint64_t B=1000000000000000000;
        std::vector<std::uint64_t> C;
        while(A>0)
        {
            auto [Q,R]=A.euclidean_division(B);
            C.push_back(static_cast<std::uint64_t>(R));
            A=std::move(Q);
        }
        std::reverse(C.begin(),C.end());
        if(C.empty()) C.push_back(0);
        for(int i=0;i<C.size();i++)
        {
            if(i) H << std::setw(18) << std::setfill('0') << C[i];
            else H << C[i];
        }
        return H;
    }


    template<std::size_t M>
    std::istream& operator>>(std::istream &H,big_int<M> &A)
    {
        std::string S;
        H >> S;
        A=0;
        for(int i=0;i<S.size();i+=18)
        {
            std::uint64_t x=0;
            for(int j=0;j<18&&i+j<S.size();j++)
                x=x*10+S[i+j]-'0';
            A=A*1000000000000000000+x;
        }
        return H;
    }
}
#endif //CPLIBRARY_BIGINT_H

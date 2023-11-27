//
// Created by ramizouari on 22/11/23.
//
#include <span>
#include "arithmetic/varint.h"
#include "arithmetic/bigint.h"

namespace cp
{
    void varint::reduce()
    {
        while(A.size() && A.back()==0)
            A.pop_back();
    }

    varint::varint(std::uint64_t x) {
        if(x) A.push_back(x);
    }

    varint varint::operator+(const varint &B) const
    {
        const auto &[P,Q]=std::minmax(A,B.A,[](auto &A,auto &B){return A.size()<B.size();});
        auto N=P.size(),M=Q.size();
        varint C;
        C.A.reserve(M+1);
        std::uint64_t carry=0;
        for(int i=0;i<N;i++)
        {
            C.A.push_back(P[i]+Q[i]);
            C.A[i]+=carry;
            carry=(C.A[i]<carry);
        }
        for(int i=N;i<M;i++)
        {
            C.A.push_back(Q[i]+carry);
            carry=(C.A[i]<carry);
        }
        if(carry) C.A.push_back(carry);
        C.reduce();
        return C;
    }

    varint &varint::operator+=(const varint &B)
    {
        if(A.size()<B.A.size()) A.resize(B.A.size());
        std::uint64_t carry=0;
        for(int i=0;i<B.A.size();i++)
        {
            A[i]+=B.A[i];
            A[i]+=carry;
            carry=(A[i]<carry);
        }
        for(int i=B.A.size();i<A.size();i++)
        {
            A[i]+=carry;
            carry=(A[i]<carry);
        }
        if(carry) A.push_back(carry);
        reduce();
        return *this;
    }

    varint varint::operator-(const varint &B) const
    {
        varint C;
        std::uint64_t carry=0;
        for(int i=0;i<std::max(A.size(),B.A.size());i++)
        {
            if(i<A.size()) C.A.push_back(A[i]);
            if(i<B.A.size()) C.A[i]-=B.A[i];
            C.A[i]-=carry;
            carry=(C.A[i]>carry);
        }
        if(carry) C.A.push_back(carry);
        C.reduce();
        return C;
    }

    varint &varint::operator-=(const varint &B)
    {
        std::uint64_t carry=0;
        for(int i=0;i<std::max(A.size(),B.A.size());i++)
        {
            if(i<A.size()) A[i]-=B.A[i];
            else A.push_back(-B.A[i]);
            A[i]-=carry;
            carry=(A[i]>carry);
        }
        if(carry) A.push_back(carry);
        reduce();
        return *this;
    }

    varint varint::operator*(const varint &B) const
    {
        varint C;
        auto N=A.size()+B.A.size();
        C.A.resize(N);
        for(int i=0;i<A.size();i++)
        {
            std::uint64_t D1=0,D2=0;
            for(int j=0;j<B.A.size();j++)
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
                C.A[i+B.A.size()] = reduce_carry(C.A[i+B.A.size()],D1);

            }
        }
        C.reduce();
        return C;
    }

    varint &varint::operator*=(const varint &B)
    {
        auto X=(*this)*B;
        A=std::move(X.A);
        return *this;
    }

    varint varint::operator<<(int k) const {
        varint C;
        C.A.resize(A.size()+k);
        for(int i=0;i<A.size();i++)
            C.A[i+k]=A[i];
        C.reduce();
        return C;
    }

    varint varint::operator>>(int k) const {
        varint C;
        C.A.resize(A.size()-k);
        for(int i=0;i<C.A.size();i++)
            C.A[i]=A[i+k];
        return C;
    }

    varint &varint::operator<<=(unsigned int k) {
        unsigned q=k/64;
        unsigned r=k&0x3F;
        A.resize(A.size()+k);
        for(int i=A.size()-1;i>=q;i--)
        {
            A[i] = A[i - q];
        }
        for(int i=0;i<k;i++)
            A[i]=0;
        reduce();
        return *this;
    }

    varint &varint::operator>>=(unsigned int k) {
        for(int i=0;i<A.size()-k;i++)
            A[i]=A[i+k];
        A.resize(A.size()-k);
        reduce();
        return *this;
    }

    std::pair<varint, varint> varint::euclidean_division(const varint &B) const
    {
        using namespace bits;
        auto N=A.size();
        auto M=B.A.size();
        std::vector<bool> X(64*N),Y(64*M);
        for(int i=0;i<N;i++) for(int j=0;j<64;j++)
                X[i*64+j]=(A[i]>>j)&1;
        for(int i=0;i<M;i++) for(int j=0;j<64;j++)
                Y[i*64+j]=(B.A[i]>>j)&1;
        std::vector<bool> Q(64*N),R;
        for(int i=64*N-1;i>=0;i--)
        {
            R<<=1;
            if(R.empty()) R.push_back(X[i]);
            else R[0]=X[i];
            while(R.size()&& !R.back()) R.pop_back();
            if(bits_compare(R,Y) >= 0)
            {
                R-=Y;
                Q[i]=true;
            }
        }
        varint C,D;
        C.A.resize((Q.size()+63)/64);
        D.A.resize((R.size()+63)/64);
        for(size_t i=0;i<Q.size();i++)
            C.A[i/64]|=static_cast<std::uint64_t>(Q[i])<<(i%64);
        for(size_t i=0;i<R.size();i++)
            D.A[i/64]|=static_cast<std::uint64_t>(R[i])<<(i%64);
        C.reduce();
        D.reduce();
        return {C,D};
    }

    varint varint::operator/(const varint &B) const
    {
        auto [X,Y]=euclidean_division(B);
        return X;
    }

    varint varint::operator%(const varint &B) const
    {
        auto [X,Y]=euclidean_division(B);
        return Y;
    }

    varint &varint::operator/=(const varint &B)
    {
        auto X=(*this)/B;
        A=std::move(X.A);
        return *this;
    }

    namespace
    {
        std::vector<std::uint64_t> karatsuba_multiplication_u64(std::span<std::uint64_t> A,std::span<std::uint64_t> B, std::uint64_t L)
        {
            if(std::min(A.size(),B.size()) <= L)
            {
                varint C,D;
                C.A.resize(A.size());
                D.A.resize(B.size());
                for(int i=0;i<A.size();i++)
                    C.A[i]=A[i];
                for(int i=0;i<B.size();i++)
                    D.A[i]=B[i];
                return (C*D).A;
            }
            else
            {
                auto N=A.size();
                auto M=B.size();
                auto K=std::max(N,M)/2;
                auto X=A.subspan(0,K),Y=A.subspan(K);
                auto U=B.subspan(0,K),V=B.subspan(K);
                auto P=karatsuba_multiplication_u64(X,U,L);
                auto Q=karatsuba_multiplication_u64(Y,V,L);
                for(int i=0;i<K;i++)
                {
                    X[i] += Y[i];
                    U[i] += V[i];
                }
                auto R=karatsuba_multiplication_u64(X,V,L);
                for(int i=0;i<K;i++)
                {
                    R[i]-=P[i];
                    R[i]-=Q[i];
                }
                std::vector<std::uint64_t> C(N+M+1);
                std::uint64_t carry=0;
                for(int i=0;i<P.size();i++)
                {
                    C[i]= reduce_carry(C[i],carry);
                    C[i]= add_with_carry(C[i],P[i],carry);
                }

                for(int i=0;i<R.size();i++)
                {
                    C[i+K]= reduce_carry(C[i+K],carry);
                    C[i+K]= add_with_carry(C[i+K],R[i]-P[i]-Q[i],carry);
                }
                for(int i=0;i<Q.size();i++)
                {
                    C[i+2*K]= reduce_carry(C[i+2*K],carry);
                    C[i+2*K]= add_with_carry(C[i+2*K],Q[i],carry);
                }
                C[N+M]=carry;
                return C;
            }
        }
    }

    varint &varint::operator%=(const varint &B)
    {
        auto X=(*this)%B;
        A=std::move(X.A);
        return *this;
    }

    std::strong_ordering varint::operator<=>(const varint &B) const
    {
        auto n=A.size(),m=B.A.size();
        while(n>0&& !A[n-1]) n--;
        while(m>0&& !B.A[m-1]) m--;
        auto p=A.size()-n;
        auto q=B.A.size()-m;
        return std::lexicographical_compare_three_way(A.rbegin()+p,A.rend(),B.A.rbegin()+q,B.A.rend());
    }

    varint varint::operator~() const     {
        varint C=*this;
        for(auto &x:C.A) x=~x;
        return C;
    }

    varint varint::operator-() const     {
        return ~(*this)+1;
    }

    varint &varint::operator++()
    {
        *this+=1;
        return *this;
    }

    varint varint::operator++(int)
    {
        varint C=*this;
        ++(*this);
        return C;
    }

    varint &varint::operator--()
    {
        *this-=1;
        return *this;
    }

    varint varint::operator--(int)
    {
        varint C=*this;
        --(*this);
        return C;
    }

    varint::operator std::uint64_t() const
    {
        return A.empty()?0:A[0];
    }

    std::ostream &operator<<(std::ostream &H, varint A)
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

    std::istream &operator>>(std::istream &H, varint &A) {
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
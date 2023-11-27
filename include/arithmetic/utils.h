//
// Created by ramizouari on 22/11/23.
//

#ifndef CPLIBRARY_UTILS_H
#define CPLIBRARY_UTILS_H
#include <cstdint>
#include <vector>

inline std::uint64_t add_with_carry(std::uint64_t x,std::uint64_t y,std::uint64_t &carry)
{
    std::uint64_t z=x+y;
    carry+=(z<x)||(z<y);
    return z;
}

inline std::uint64_t reduce_carry(std::uint64_t x,std::uint64_t &carry)
{
    std::uint64_t z=x+carry;
    carry=(z<x);
    return z;
}

inline std::uint64_t add_with_carry(std::uint64_t x,std::uint64_t y, std::uint64_t z,std::uint64_t &carry)
{
    return add_with_carry(add_with_carry(x,y,carry),z,carry);
}

namespace bits
{
    template<std::size_t N>
    std::strong_ordering operator<=>(const std::bitset<N> &A,const std::bitset<N> &B)
    {
        for(int i=N-1;i>=0;i--)
        {
            if(A[i]<B[i]) return std::strong_ordering::less;
            if(A[i]>B[i]) return std::strong_ordering::greater;
        }
        return std::strong_ordering::equal;
    }

    template<std::size_t N>
    std::bitset<N> operator+(const std::bitset<N> &A,const std::bitset<N> &B)
    {
        std::bitset<N> C;
        std::uint64_t carry=0;
        for(int i=0;i<N;i++)
        {
            C[i]=A[i]^B[i]^carry;
            carry=(A[i]&B[i])||(A[i]&carry)||(B[i]&carry);
        }
        return C;
    }

    inline std::vector<bool> operator+(const std::vector<bool> &A,const std::vector<bool> &B)
    {
        auto [P,Q]=std::minmax(A,B,[](auto &A,auto &B){return A.size()<B.size();});
        auto n=P.size(),m=Q.size();
        std::vector<bool> C(std::max(n,m));
        bool carry=false;
        for(int i=0;i<n;i++)
        {
            C[i]=P[i]^Q[i]^carry;
            carry=(P[i]&&Q[i])||(P[i]&&carry)||(Q[i]&&carry);
        }
        for(int i=n;i<m;i++)
        {
            C[i]=Q[i]^carry;
            carry=Q[i]&&carry;
        }
        if(carry) C.push_back(carry);
        return C;
    }

    inline std::vector<bool>& operator+=(std::vector<bool> &A,const std::vector<bool> &B)
    {
        return A=A+B;
    }

    inline std::vector<bool> operator~(const std::vector<bool> &A)
    {
        std::vector<bool> B;
        for(auto x:A)
            B.push_back(!x);
        return B;
    }

    inline std::vector<bool> operator-(const std::vector<bool> &A,const std::vector<bool> &B)
    {
        auto r=std::max(A.size(),B.size());
        std::vector<bool> C(A),D(~B);
        C.resize(r,false);
        D.resize(r,true);
        D+=std::vector<bool>(1,true);
        auto Z= C+D;
        Z.resize(r);
        while(Z.size()&& !Z.back())
            Z.pop_back();
        return Z;
    }


    inline std::vector<bool>& operator-=(std::vector<bool> &A,const std::vector<bool> &B)
    {
        return A=A-B;
    }

    inline std::vector<bool> operator<<(const std::vector<bool> &A,size_t k)
    {
        std::vector<bool> B;
        for(int i=0;i<k;i++)
            B.push_back(false);
        for(auto x:A)
            B.push_back(x);
        return B;
    }

    inline std::vector<bool> operator>>(const std::vector<bool> &A,size_t k)
    {
        std::vector<bool> B;
        for(int i=0;i<A.size()-k;i++)
            B.push_back(A[i+k]);
        return B;
    }

    inline std::vector<bool>& operator<<=(std::vector<bool> &A,size_t k)
    {
        A.resize(A.size()+k);
        for(int i=A.size()-1;i>=k;i--)
            A[i]=A[i-k];
        for(int i=0;i<k;i++)
            A[i]=false;
        return A;
    }

    inline std::vector<bool>& operator>>=(std::vector<bool> &A,size_t k)
    {
        for(int i=0;i<A.size()-k;i++)
            A[i]=A[i+k];
        A.resize(A.size()-k);
        return A;
    }

    template<std::size_t N>
    std::bitset<N> operator-(const std::bitset<N> &A,const std::bitset<N> &B)
    {
        return A+~B+1;
    }

    template<std::size_t N>
    std::bitset<N>& operator+=(std::bitset<N> &A,const std::bitset<N> &B)
    {
        bool carry=false;
        for(int i=0;i<N;i++)
        {
            auto tmp=A[i]^B[i]^carry;
            carry=(A[i]&&B[i])||(A[i]&&carry)||(B[i]&&carry);
            A[i]=tmp;
        }
        return A;
    }

    template<std::size_t N>
    std::bitset<N>& operator-=(std::bitset<N> &A,const std::bitset<N> &B)
    {
        auto C=~B+std::bitset<N>(1);
        return A+=C;
    }

    inline std::weak_ordering bits_compare(const std::vector<bool> &A,const std::vector<bool> &B)
    {
        auto i=A.size();
        auto j=B.size();
        while(i>0&& !A[i-1]) i--;
        while(j>0&& !B[j-1]) j--;
        auto p=A.size()-i;
        auto q=B.size()-j;
        if(i!=j)
            return i<=>j;
        return std::lexicographical_compare_three_way(A.rbegin()+p,A.rend(),B.rbegin()+q,B.rend());
    }
}

#endif //CPLIBRARY_UTILS_H

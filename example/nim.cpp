#include <iostream>
#include <vector>
#include <chrono>
#include "polynomial/polynomial.h"
#include "signals/ntt.h"
#include "polynomial/ring_extension.h"
#include "polynomial/fast_polynomial.h"
#include "polynomial/special_polynomials.h"
#include "polynomial/fast_general_polynomial_multiplication.h"
#include "arithmetic/bigint.h"
#include "arithmetic/fast_arithmetic.h"

namespace cp::signals
{

    template<typename R>
    struct cyclic_extension
    {
        std::vector<R> p;
    public:
        cyclic_extension(const std::vector<R> & _p):p(_p)
        {
        }
        cyclic_extension(R _p):p(1,_p)
        {
        }
        cyclic_extension()=default;

        cyclic_extension& operator+=(const cyclic_extension &O)
        {
            int r=std::max(p.size(),O.p.size());
            p.resize(r);
            for(int i=0;i<O.p.size();i++)
                p[i]+=O.p[i];
            return *this;
        }

        cyclic_extension& operator-=(const cyclic_extension &O)
        {
            int r=std::max(p.size(),O.p.size());
            p.resize(r);
            for(int i=0;i<O.p.size();i++)
                p[i]-=O.p[i];
            return *this;
        }

        cyclic_extension operator*(const cyclic_extension &O) const
        {
            int n=p.size(),m=O.p.size();
            cyclic_extension q;
            auto r=std::max(n,m);
            q.p.resize(r);
            for(int i=0;i<n;i++) for(int j=0;j<m;j++)
                    q.p[(i+j)%r]+=p[i]*O.p[j];
            if(std::has_single_bit<unsigned>(n)) for(int i=r/2;i<r;i++)
            {
                q.p[i - r / 2] -= q.p[i];
                q.p[i] = 0;
            }

            else
                throw std::runtime_error("Not implemented");

            return q;
        }

        int degree() const
        {
            return p.size();
        }

        cyclic_extension operator*(const R &O) const
        {
            cyclic_extension q;
            q.p.resize(degree());
            for(int i=0;i<degree();i++)
                q.p[i]=p[i]*O;
            return q;
        }

        cyclic_extension operator+(const cyclic_extension &O) const
        {
            cyclic_extension q;
            q.p.resize(std::max(p.size(),O.p.size()));
            for(int i=0;i<p.size();i++)
                q.p[i]=p[i];
            for(int i=0;i<O.p.size();i++)
                q.p[i]+=O.p[i];
            return q;
        }

        cyclic_extension operator-(const cyclic_extension &O) const
        {
            cyclic_extension q;
            q.p.resize(std::max(p.size(),O.p.size()));
            for(int i=0;i<p.size();i++)
                q.p[i]=p[i];
            for(int i=0;i<O.p.size();i++)
                q.p[i]-=O.p[i];
            return q;
        }

        cyclic_extension operator-() const
        {
            cyclic_extension q;
            q.p.resize(p.size());
            for(int i=0;i<p.size();i++)
                q.p[i]=-p[i];
            return q;
        }

        cyclic_extension operator*=(const cyclic_extension &O)
        {
            auto r=(*this)*O;
            p.swap(r.p);
            return *this;
        }

        cyclic_extension operator*=(const R &O)
        {
            for(auto &s:p)
                s*=O;
            return *this;
        }

        cyclic_extension operator+=(const R &O)
        {
            p[0]+=O;
            return *this;
        }

        cyclic_extension operator-=(const R &O)
        {
            p[0]-=O;
            return *this;
        }

        cyclic_extension operator+(const R &O) const
        {
            cyclic_extension q;
            q.p.resize(p.size());
            for(int i=0;i<p.size();i++)
                q.p[i]=p[i];
            q.p[0]+=O;
            return q;
        }

        cyclic_extension operator-(const R &O) const
        {
            cyclic_extension q;
            q.p.resize(p.size());
            for(int i=0;i<p.size();i++)
                q.p[i]=p[i];
            q.p[0]-=O;
            return q;
        }

    };

    template<typename R>
    cyclic_extension<R> operator*(const R &a, const cyclic_extension<R> &b)
    {
        return b*a;
    }

    template<typename R>
    cyclic_extension<R> operator+(const R &a, const cyclic_extension<R> &b)
    {
        return b+a;
    }

    template<typename R>
    cyclic_extension<R> operator-(const R &a, const cyclic_extension<R> &b)
    {
        return -b+a;
    }

    template<typename R>
    void inplace_ring_ntt2(cp::linalg::tensor_view<cyclic_extension<R>,1> & a, bool inverse=false, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        int n = a.size();
        for (int i = 1, j = 0; i < n; i++)
        {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                std::swap(a(i), a(j));
        }
        std::vector<cyclic_extension<R>> W;
        auto deg=a(0).degree();
        std::vector<R> p(deg);
        p[inverse?deg-1:1]=1;
        W.push_back(cyclic_extension<R>(p));
        for (int len = 2; len <= n; len <<= 1)
            W.push_back(W.back()*W.back());
        W.pop_back();
        for (int len = 2,r=W.size()-1; len <= n; len <<= 1,r--)
        {
            auto wlen = W[r];
            for (int i = 0; i < n; i += len)
            {
                cyclic_extension<R> w(1);
                for (int j = 0; j < len / 2; j++)
                {
                    cyclic_extension<R> u = a(i+j), v = a(i+j+len/2) * w;
                    a(i+j) = u + v;
                    a(i+j+len/2) = u - v;
                    w *= wlen;
                }
            }
        }
    }

    template<typename R>
    void inplace_ring_ntt2(cp::linalg::tensor_view<cyclic_extension<R>,1> && a, bool inverse=false, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        inplace_ring_ntt2(a,inverse,normalized);
    }
}

constexpr cp::integer M=998244353;
using IK=cp::cyclic<M>;

cp::varint factorial(int n)
{
    cp::varint P=1;
    for(int i=1;i<=n;i++)
        P*=i;
    return P;
}

cp::varint fast_factorial(int n)
{
    std::priority_queue<cp::varint,std::vector<cp::varint>,std::greater<>> Q;
    for(int i=1;i<=n;i++)
        Q.push(i);
    while(Q.size()>1)
    {
        auto a=Q.top();
        Q.pop();
        auto b=Q.top();
        Q.pop();
        if(std::min(a.A.size(), b.A.size()) < 1)
            Q.push(a*b);
        else
            Q.push(cp::fast_multiplication(a,b));
    }
    return Q.top();
}

template<typename R,int r>
void print(std::ostream& H,const std::vector<cp::details::cyclic_ring<R,r>> &A)
{
    for(auto &s:A)
    {
        std::cout << "[ ";
        for(auto &t:s.p)
            std::cout << t << ", ";
        std::cout << "]\n";
    }

}
template<typename R,int r,bool inverse>
void print(std::ostream& H,const cp::details::cyclic_ring_extension<R,r,inverse>&s)
{
    H << "[ ";
    for(auto a:s.data())
    {
        H << "( ";
        for(auto b:a.p)
            H << b << ", ";
        H << "); ";
    }
    H << "]\n";

}


template<typename R,int r,bool inverse>
void print(std::ostream& H,const std::vector<cp::details::cyclic_ring_extension<R,r,inverse>> &A)
{
    for(auto &s:A)
    {
        H << "[ ";
        for(auto a:s.data())
        {
            H << "( ";
            for(auto b:a.p)
                H << b << ", ";
            H << "); ";
        }
        H << "]\n";
    }

}

using R=std::uint64_t;

int main()
{
    auto F= std::make_shared<cp::factoriser>(1e6);
    int T;
    std::cin >> T;
    while(T--)
    {
        int n;
        constexpr int s=3;
        cp::default_factoriser_t::default_factoriser=F;
        std::cin >> n;
        std::vector<R> T(n,1);
        auto Op=[](const cp::polynomial<R> &a,const  cp::polynomial<R> &b) -> cp::polynomial<R>{return cp::karatsuba_multiplication(a,b);};
        std::function Fn=Op;
        auto t1=std::chrono::high_resolution_clock::now();
        int block=128;
#ifndef NDEBUG
        block=4;
#endif
        auto [Y,m]=cp::fast_general_polynomial_multiplication_p2(T,T,block);
        auto t2=std::chrono::high_resolution_clock::now();
        auto fn=[](const auto &a,const  auto &b) {return cp::karatsuba_multiplication(a,b);};

        auto Z= cp::karatsuba_multiplication(cp::polynomial(T),cp::polynomial(T));
        auto t3=std::chrono::high_resolution_clock::now();
        std::cout << "Multiplier: " << m << "\n";
        for(int i=0;i<Y.size();i++)
            std::cout << Y[i]/(1<<m) << " ";
        std::cout << "\n";
        std::cout << "\n";
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "\n";
        std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t3-t2).count() << "\n";
    }

}
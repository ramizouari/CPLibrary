//
// Created by ASUS on 01/12/2021.
//

#ifndef ACPC_PREPARATION_FAST_POLY_H
#define ACPC_PREPARATION_FAST_POLY_H
#include "fft.h"
#include "ring_extension.h"
#include "sparse_polynomial.h"
#include "algebra/binary_operation.h"
#include "data_structures/range_queries.h"
#include "data_structures/fixed/segment_tree.h"
#include "signals/fft.h"
#include "linear_algebra/view.h"
#include "data_structures/dynamic/segment_tree.h"

namespace cp
{

    /**
 * @brief Karatsuba multiplication
* @details Applies Karatsuba multiplication between two polynomials
* @Requirements
* None
*/

    template<typename R>
    polynomial<R> karatsuba_multiplication(const polynomial<R> &p,const polynomial<R> &q)
    {
        constexpr int L=75;
        if(std::min(p.degree(),q.degree())<=L)
            return p*q;
        polynomial<R> a1,b1,a2,b2;
        int n=p.degree(),m=q.degree(),r=std::max(n,m)+1;
        std::vector<R> &u1=a1.data(),&u2=a2.data(),
                &v1=b1.data(),&v2=b2.data();
        u1.resize(std::min(n+1,r/2));
        u2.resize(std::min(m+1,r/2));
        v1.resize(std::max(n+1-r/2,0));
        v2.resize(std::max(m+1-r/2,0));
        for(int i=0;i<u1.size();i++)
            u1[i]=p[i];
        for(int i=0;i<u2.size();i++)
            u2[i]=q[i];
        for(int i=0;i<v1.size();i++)
            v1[i]=p[i+r/2];
        for(int i=0;i<v2.size();i++)
            v2[i]=q[i+r/2];
        polynomial<R> r1= karatsuba_multiplication(a1,a2),
                r3= karatsuba_multiplication(b1,b2),
                t=karatsuba_multiplication(a1+b1,a2+b2),
                r2=t-r1-r3;
        polynomial<R> h;
        int s=r-r%2;
        auto &c=h.data();
        c.resize(n+m+1);
        for(int i=0;i<=r1.degree();i++)
            c[i]+=r1[i];
        for(int i=0;i<=r2.degree();i++)
            c[i+r/2]+=r2[i];
        for(int i=0;i<=r3.degree();i++)
            c[i+s]+=r3[i];
        return h;
    }

    template<typename R>
    sparse_polynomial<R> karatsuba_multiplication(const sparse_polynomial<R>& p, const sparse_polynomial<R>& q)
    {
        constexpr int recursion_limit = 30;
        if (std::min(p.size(), q.size()) <= recursion_limit)
            return p * q;
        sparse_polynomial<R> a1, b1, a2, b2;
        int n = p.degree(), m = q.degree(), r = std::max(n, m) + 1;
        const auto &mapper1 = static_cast<const std::map<int, R>&>(p),&mapper2=static_cast<const std::map<int, R>&>(q);
        auto it1 = mapper1.begin(),it2=mapper2.begin();
        for (; it1!=mapper1.end()  && it1->first<  r / 2; ++it1)
            a1[it1->first] = it1->second;
        for (; it2 != mapper2.end() && it2->first < r / 2; ++it2)
            a2[it2->first] = it2->second;
        for (; it1 != mapper1.end(); ++it1)
            b1[it1->first-r/2] = it1->second;
        for (; it2 != mapper2.end(); ++it2)
            b2[it2->first-r/2] = it2->second;
        sparse_polynomial<R> r1 = karatsuba_multiplication(a1, a2),
                r3 = karatsuba_multiplication(b1, b2),
                t = karatsuba_multiplication(a1 + b1, a2 + b2),
                r2 = t - r1 - r3;
        sparse_polynomial<R> h;
        int s = r - r % 2;
        auto& c = static_cast<std::map<int,R>&>(h);
        c = r1;
        for (auto [k, w] : static_cast<std::map<int, R>&>(r2))
            c[k + r / 2] += w;
        for (auto [k, w] : static_cast<std::map<int, R>&>(r3))
            c[k + s] += w;
        return h;
    }

    template<typename R>
    struct polynomial_operation;
    template<typename R>
    struct polynomial_operation<polynomial<R>> : public binary_operation<polynomial<R>>
    {
        std::shared_ptr<binary_operation<std::vector<R>>> m_op;
        std::shared_ptr<invertible_operation<R>> m_inv;
        polynomial_operation(std::shared_ptr<binary_operation<std::vector<R>>> _op,std::shared_ptr<invertible_operation<R>> _inv):m_op(_op),m_inv(_inv){}
        polynomial_operation(std::shared_ptr<binary_operation<std::vector<R>>> _op):m_op(_op),m_inv(nullptr){}
        polynomial<R> reduce(const polynomial<R>&a,const polynomial<R>&b) const override
        {
            return m_op->reduce(a.data(),b.data());
        }

        R inv(const R&a) const
        {
            return m_inv->inv(a);
        }

        polynomial<R> neutral_element() const override
        {
            return polynomial<R>(m_op->neutral_element());
        }
        std::shared_ptr<binary_operation<std::vector<R>>> underlying_operator() const
        {
            return m_op;
        }
        std::shared_ptr<invertible_operation<R>> scalar_inverter() const
        {
            return m_inv;
        }
    };

    template<typename R>
    struct karatsuba_multiplies_t : public polynomial_operation<R>
    {
        R reduce(const R&a,const R&b) const override
        {
            return karatsuba_multiplication(a,b);
        }
    };

    template<typename R>
    struct fast_multiplies_t;

    template<typename R>
    struct fast_multiplies_t<std::vector<R>> : public binary_operation<std::vector<R>>
    {
        signals::radix2_fft<R> fft;
        std::vector<R> reduce(const std::vector<R>&a,const std::vector<R>&b) const override
        {
            if(a.size()==0 || b.size()==0)
                return {};
            std::vector<R> A(a),B(b);
            auto r=std::bit_ceil<unsigned>(A.size()+B.size()-1);
            A.resize(r);
            B.resize(r);
            linalg::vector_view U(A),V(B);
            fft.transform(U,false,signals::FFTNormalization::None);
            fft.transform(V,false,signals::FFTNormalization::None);
            for(int i=0;i<r;i++)
                A[i]*=B[i];
            fft.transform(U,true,signals::FFTNormalization::Normalized);
            A.resize(a.size()+b.size()-1);
            return A;
        }
        inline static std::vector<R> neutral=std::vector<R>{1};
        std::vector<R> neutral_element() const override
        {
            return neutral;
        }
    };

    template<std::floating_point R>
    struct fast_multiplies_t<std::vector<R>> : public binary_operation<std::vector<R>>
    {
        signals::radix2_fft<std::complex<R>> fft;
        std::vector<R> reduce(const std::vector<R>&a,const std::vector<R>&b) const override
        {
            if(a.size()==0 || b.size()==0)
                return {};
            auto r=std::bit_ceil<unsigned>(a.size()+b.size()-1);
            std::vector<std::complex<R>> A(r);
            for(int i=0;i<a.size();i++)
            {
                A[i].real(a[i]);
                A[i].imag(b[i]);
            }
            linalg::vector_view U(A);
            fft.transform(U,false,signals::FFTNormalization::None);
            std::vector<std::complex<R>> C(r);
            for(unsigned i=0;i<r;i++)
            {
                auto j=(r-i)&(r-1);
                C[i]=A[i]*A[i]-std::conj(A[j]*A[j]);
                C[i]*=std::complex<R>(0,-0.25);
            }
            linalg::vector_view V(C);
            fft.transform(V,true,signals::FFTNormalization::Normalized);
            std::vector<R> C_real(r);
            for(int i=0;i<r;i++)
                C_real[i]=C[i].real();
            C_real.resize(a.size()+b.size()-1);
            return C_real;
        }
        inline static std::vector<R> neutral=std::vector<R>{1};
        std::vector<R> neutral_element() const override
        {
            return neutral;
        }
    };

    template<typename R>
    struct fast_multiplies_t<polynomial<R>> : public binary_operation<polynomial<R>>
    {
        fast_multiplies_t<std::vector<R>> fft;
        polynomial<R> reduce(const polynomial<R>&a,const polynomial<R>&b) const override
        {
            return polynomial<R>(fft.reduce(a.data(),b.data()));
        }
        inline static polynomial<R> neutral=polynomial<R>{1};
    };

    template<typename R>
    std::vector<R> fast_multiplication(const std::vector<R> &A,const std::vector<R> &B)
    {
        static fast_multiplies_t<std::vector<R>> multiplies;
        return multiplies.reduce(A,B);
    }

    template<typename R>
    polynomial<R> fast_multiplication(const polynomial<R> &A,const polynomial<R> &B)
    {
        static fast_multiplies_t<polynomial<R>> multiplies;
        return multiplies.reduce(A,B);
    }

    template<typename Real=real, cp::integer m>
    std::vector<cp::cyclic<m>> fast_modular_multiplication_real(const std::vector<cp::cyclic<m>> &a,const std::vector<cp::cyclic<m>> &b)
    {
        if(a.size()==0)
            return {};
        using namespace cp;
        std::array<std::vector<Real>,2> A,B;
        integer block=std::ceil(std::sqrt(a.front().modulus()));
        for(int i=0;i<a.size();i++)
        {
            auto [q,r]=std::div(static_cast<integer>(a[i]),block);
            A[0].push_back(r);
            A[1].push_back(q);
        }
        for(int i=0;i<b.size();i++)
        {
            auto [q,r]=std::div(static_cast<integer>(b[i]),block);
            B[0].push_back(r);
            B[1].push_back(q);
        }
        std::array<std::vector<Real>,4> C;
        C[0]=fast_multiplication(A[0],B[0]);
        C[1]=fast_multiplication(A[0],B[1]);
        C[2]=fast_multiplication(A[1],B[0]);
        C[3]=fast_multiplication(A[1],B[1]);
        std::vector<cyclic<m>> R(a.size()+b.size()-1);
        for(int i=0;i<R.size();i++)
        {
            integer x=0;
            x+=std::llround(C[0][i]);
            x+=std::llround(C[1][i])%m*block;
            x+=std::llround(C[2][i])%m*block%m;
            x+=std::llround(C[3][i])%m*(block*block)%m;
            R[i]=x;
        }
        return R;
    }


    template<typename R>
    std::vector<R> formal_inv_2(const std::vector<R> &A,int m)
    {
        if(m==1)
            return {R(1)/A.front()};
        auto B=A;
        for(int i=1;i<A.size();i+=2)
            B[i]=-B[i];
        auto C= fast_multiplication(A,B);
        std::vector<R> T;
        T.resize(m/2);
        for(int i=0;i<T.size() && 2*i < C.size();i++)
            T[i]=C[2*i];
        auto S=formal_inv_2(T,m/2);
        std::vector<R> Q;
        Q.resize(m);
        for(int i=0;i<m/2;i++)
            Q[2*i]=S[i];
        return fast_multiplication(B, Q);
    }

    template<typename R>
    std::vector<R> formal_inv_2(const std::vector<R> &A,int m,std::shared_ptr<binary_operation<std::vector<R>>> multiplies,std::shared_ptr<invertible_operation<R>> inv)
    {
        if(m==1)
            return {inv->inv(A.front())};
        auto B=A;
        for(int i=1;i<A.size();i+=2)
            B[i]=-B[i];
        auto C= multiplies->reduce(A,B);
        std::vector<R> T;
        T.resize(m/2);
        for(int i=0;i<T.size() && 2*i < C.size();i++)
            T[i]=C[2*i];
        auto S=formal_inv_2(T,m/2,multiplies,inv);
        std::vector<R> Q;
        Q.resize(m);
        for(int i=0;i<m/2;i++)
            Q[2*i]=S[i];
        return multiplies->reduce(B, Q);
    }

    template<typename R>
    polynomial<R> formal_inv_2(const polynomial<R> &A,int m)
    {
        return formal_inv_2(A.data(),m);
    }

    template<typename R>
    polynomial<R> formal_inv_2(const polynomial<R> &A,int m,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        return formal_inv_2(A.data(),m,multiplies->underlying_operator(),multiplies->scalar_inverter());
    }

    template<typename R>
    std::vector<R> formal_inv(const std::vector<R> &A,int m)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m));
        C.resize(m);
        return C;
    }

    template<typename R>
    std::vector<R> formal_inv(const std::vector<R> &A,int m,std::shared_ptr<binary_operation<std::vector<R>>> multiplies,std::shared_ptr<invertible_operation<R>> inv)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m),multiplies,inv);
        C.resize(m);
        return C;
    }

    template<typename R>
    polynomial<R> formal_inv(const polynomial<R> &A,int m)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m));
        C.data().resize(m);
        return C;
    }

    template<typename R>
    polynomial<R> formal_inv(const polynomial<R> &A,int m,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        auto C=formal_inv_2(A,std::bit_ceil<unsigned>(m),multiplies, multiplies->scalar_inverter());
        C.data().resize(m);
        return C;
    }

    template<typename R,int nilpotence>
    nilpotent_extension<R,nilpotence> fast_inv(const nilpotent_extension<R,nilpotence> &A)
    {
        return nilpotent_extension<R,nilpotence>(formal_inv(A.p,nilpotence));
    }

    template<typename R>
    d_nilpotent_extension<R> fast_inv(const d_nilpotent_extension<R> &A)
    {
        return d_nilpotent_extension<R>(formal_inv(A.p),nilpotence_t{A.nilpotence});
    }


    template<typename R>
    std::vector<R> fast_division(std::vector<R> A,std::vector<R> Q)
    {
        if(A.size()<Q.size())
            return {};
        int m=A.size()-Q.size()+1;
        std::reverse(A.begin(),A.end());
        std::reverse(Q.begin(),Q.end());
        auto P= fast_multiplication(A, formal_inv(Q,m));
        P.resize(m);
        std::reverse(P.begin(),P.end());
        return P;
    }

    template<typename R>
    std::vector<R> fast_division(std::vector<R> A,std::vector<R> Q,std::shared_ptr<binary_operation<std::vector<R>>> multiplies,std::shared_ptr<invertible_operation<R>> inv)
    {
        if(A.size()<Q.size())
            return {};
        int m=A.size()-Q.size()+1;
        std::reverse(A.begin(),A.end());
        std::reverse(Q.begin(),Q.end());
        auto P= multiplies->reduce(A, formal_inv(Q,m,multiplies,inv));
        P.resize(m);
        std::reverse(P.begin(),P.end());
        return P;
    }

    template<typename R>
    polynomial<R> fast_division(const polynomial<R> &A,const polynomial<R> &B)
    {
        return fast_division(A.data(),B.data());
    }

    template<typename R>
    polynomial<R> fast_division(const polynomial<R> &A,const polynomial<R> &B, std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        return fast_division(A.data(),B.data(), multiplies->underlying_operator(),multiplies->scalar_inverter());
    }


    template<typename R>
    polynomial<R> fast_mod(const polynomial<R>&A,const polynomial<R>& B)
    {
        auto P= fast_division(A,B);
        auto Z= A - fast_multiplication(B,P);
        Z.data().resize(B.degree());
        return Z;
    }

    template<typename R>
    polynomial<R> fast_mod(const polynomial<R>&A, const polynomial<R>& B, std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        auto P= fast_division(A,B,multiplies);
        auto Z= A - multiplies->reduce(B,P);
        Z.data().resize(B.degree());
        return Z;
    }

    template<typename R>
    std::pair<polynomial<R>,polynomial<R>> fast_euclidean_division(const polynomial<R> &A,const polynomial<R>& B)
    {
        auto P= fast_division(A,B);
        auto Q=A- fast_multiplication(B,P);
        Q.data().resize(B.degree());
        return std::make_pair(P,Q);
    }

    template<typename R>
    polynomial<R> fast_gcd(const polynomial<R> &A,const polynomial<R> &B)
    {
        if(B==R{})
            return A;
        return fast_gcd(B,fast_mod(A,B).reduce());
    }

    template<typename R>
    polynomial<R> fast_polynomial_expansion(const std::vector<R> &X)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::fixed::segment_tree<fast_multiplies_t<polynomial<R>>> S(P);
        return S.S[0][0];
    }

    template<typename R>
    polynomial<R> fast_polynomial_expansion(const std::vector<R> &X,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::dynamic::segment_tree S(P,multiplies);
        return S.S[0][0];
    }


    template<typename R>
    std::vector<R> fast_multi_evaluation(const polynomial<R> &A,const std::vector<R> &X)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::fixed::segment_tree<fast_multiplies_t<polynomial<R>>> S(P);
        std::vector<polynomial<R>> Z(1<<(S.h-1));
        Z[0]=fast_mod(A,S.S[0][0]);
        for(int i=1;i<S.h;i++)
            for(int j=(1<<i)-1;j>=0;j--)
                Z[j]=fast_mod(Z[j>>1],S.S[i][j]);
        std::vector<R> Y;
        Y.reserve(n);
        for(int i=0;i<n;i++)
            Y.push_back(Z[i](R{}));
        return Y;
    }

    template<typename R>
    std::vector<R> fast_multi_evaluation(const polynomial<R> &A,const std::vector<R> &X,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::dynamic::segment_tree S(P,std::dynamic_pointer_cast<binary_operation<polynomial<R>>>(multiplies));
        std::vector<polynomial<R>> Z(1<<(S.h-1));
        Z[0]=fast_mod(A,S.S[0][0],multiplies);
        for(int i=1;i<S.h;i++)
            for(int j=(1<<i)-1;j>=0;j--)
                Z[j]=fast_mod(Z[j>>1],S.S[i][j],multiplies);
        std::vector<R> Y;
        Y.reserve(n);
        for(int i=0;i<n;i++)
            Y.push_back(Z[i](R{}));
        return Y;
    }

    template<typename R>
    polynomial<R> fast_interpolation(const std::vector<R> &X,const std::vector<R> &Y)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::fixed::segment_tree<fast_multiplies_t<polynomial<R>>> S(P);
        std::vector<std::vector<polynomial<R>>> Z(S.h);
        for(int i=0;i<S.h;i++)
            Z[i].resize(1<<i);
        Z[0][0]=polynomial<R>(Y);
        for(int i=1;i<S.h;i++)
            for(int j=0;j<Z[i].size();j++)
                Z[i][j]=fast_mod(Z[i-1][j>>1],S.S[i][j]);
        return Z.back()[0];
    }

    template<typename R>
    polynomial<R> fast_interpolation(const std::vector<R> &X,const std::vector<R> &Y,std::shared_ptr<polynomial_operation<polynomial<R>>> multiplies)
    {
        int n=X.size();
        std::vector<polynomial<R>> P(X.size());
        for(int i=0;i<n;i++)
            P[i]=polynomial<R>({-X[i],1});
        data_structures::dynamic::segment_tree S(P,std::dynamic_pointer_cast<binary_operation<polynomial<R>>>(multiplies));
        auto P_=S.S[0][0].derivative();
        auto Y_=fast_multi_evaluation(P_,X,multiplies);
        std::vector<polynomial<R>> Z(1<<(S.h-1));
        for(int i=0;i<Y_.size();i++)
            Z[i]=multiplies->reduce(polynomial<R>(Y[i]),polynomial<R>(multiplies->inv(Y_[i])));
        for(int i=S.h-2;i>=0;i--) for(int j=0;j<(1<<i);j++)
            Z[j]=multiplies->reduce(Z[2*j],S.S[i+1][2*j+1])+multiplies->reduce(Z[2*j+1],S.S[i+1][2*j]);
        return Z[0];
    }

    template<typename R,integer m>
    std::vector<cyclic<m>> fast_complex_modular_multiplication(const std::vector<cyclic<m>> &A,const std::vector<cyclic<m>> &B)
    {
        signals::radix2_fft<std::complex<R>> fft;
        if(A.empty())
            return {};
        auto n=std::bit_ceil(A.size()+B.size()-1);
        integer block=std::sqrt(A.front().modulus());
        std::vector<cyclic<m>> P(n),Q(n);
        for(auto a:A)
        {
            auto [x,y]=std::div(static_cast<integer>(a),block);
            P[x]=std::complex<R>(x,y);
        }
        for(auto b:B)
        {
            auto [x,y]=std::div(static_cast<integer>(b),block);
            Q[x]=std::complex<R>(x,y);
        }
        linalg::vector_view U(P),V(Q);
        fft.transform(U,false,signals::FFTNormalization::None);
        fft.transform(V,false,signals::FFTNormalization::None);

    }
}
#endif //ACPC_PREPARATION_FFT_H

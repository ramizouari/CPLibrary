//
// Created by ramizouari on 11/12/23.
//

#ifndef CPLIBRARY_FAST_GENERAL_POLYNOMIAL_MULTIPLICATION_H
#define CPLIBRARY_FAST_GENERAL_POLYNOMIAL_MULTIPLICATION_H
#include "special_polynomials.h"
#include "signals/fft.h"

/**
 * @brief This file contains an implementation of the fast general polynomial multiplication algorithm
 * @details This algorithm is a polynomial multiplication algorithm that works for any ring. It runs in O( n log(n) log(log(n)) ) time
 * */

namespace cp
{
    namespace details
    {

        template<std::signed_integral A,std::integral B=A>
        typename std::common_type<A,B>::type mod_operator(A a,B b)
        {
            return (a%b+b)%b;
        }

        template<std::integral I>
        struct circular_ratio
        {
            I p,q;
            circular_ratio(I _p=0,I _q=0):p(_p),q(_q)
            {
            }

            circular_ratio operator+(const circular_ratio &O) const
            {
                auto r=std::lcm(q,O.q);
                if(!r) return circular_ratio(p+O.p);
                return circular_ratio(mod_operator(p*(r/q)+O.p*(r/O.q),r),r);
            }

            circular_ratio operator-(const circular_ratio &O) const
            {
                auto r=std::lcm(q,O.q);
                if(!r) return circular_ratio(p-O.p);
                return circular_ratio(mod_operator(p*(r/q)-O.p*(r/O.q),r),r);
            }

            circular_ratio& operator+=(const circular_ratio &O)
            {
                auto r=std::lcm(q,O.q);
                if(!r) p+=O.p;
                else
                {
                    p*=r/q;
                    p+=O.p*(r/O.q);
                    p=mod_operator(p,r);
                    q=r;
                }
                return *this;
            }

            circular_ratio& operator-=(const circular_ratio &O)
            {
                auto r=std::lcm(q,O.q);
                if(!r) p-=O.p;
                else
                {
                    p*=r/q;
                    p-=O.p*(r/O.q);
                    p=mod_operator(p,r);
                    q=r;
                }
                return *this;
            }

            circular_ratio& operator*=(const I& O)
            {
                p*=O;
                if(q)
                    p=mod_operator(p,q);
                return *this;
            }

            circular_ratio operator*(const I& O) const
            {
                auto r=*this;
                return r*=O;
            }

            circular_ratio operator-() const
            {
                if(!q)
                    return circular_ratio(-p,q);
                else return circular_ratio(mod_operator(-p,q),q);
            }
        };

        struct root_of_unity
        {
            circular_ratio<integer> xi;
            explicit root_of_unity(integer _o=0,integer m=0):xi(_o,m)
            {
            }
            root_of_unity(const circular_ratio<integer > &_xi):xi(_xi)
            {
            }

            root_of_unity operator*(const root_of_unity &O) const
            {
                return xi+O.xi;
            }

            root_of_unity& operator*=(const root_of_unity &O)
            {
                xi+=O.xi;
                return *this;
            }

            root_of_unity operator/(const root_of_unity &O) const
            {
                return xi-O.xi;
            }

            root_of_unity& operator/=(const root_of_unity &O)
            {
                xi-=O.xi;
                return *this;
            }

            root_of_unity operator-() const
            {
                return -xi;
            }

            root_of_unity& operator^=(integer n)
            {
                xi*=n;
                return *this;
            }

            root_of_unity operator^(integer n) const
            {
                auto r=*this;
                return r^=n;
            }

            [[nodiscard]] integer position(integer order) const
            {
                if(!xi.q)
                    throw std::invalid_argument("Not a root of unity");
                auto [p,r] =std::div(xi.p*order,xi.q);
                if(r)
                    throw std::invalid_argument("Not a root of unity");
                return p;
            }
        };

        template<typename R,int s>
        struct cyclic_ring
        {
            inline static factoriser F{s+1};
            inline static polynomial<R> phi_s{},Z_s,H1{},H2{};
            inline static std::vector<std::size_t> indexes;
            std::array<R,s> p;
            cyclic_ring(const polynomial<R> &P):p{}
            {
                for(int i=0;i<std::min<int>(s,P.data().size());i++)
                    p[i]=P[i];
            }
            cyclic_ring(const R &a=0):p{}
            {
                p[0]=a;
            }

            template<std::integral A>
            cyclic_ring(const A &a):p{}
            {
                p[0]=a;
            }

            cyclic_ring& operator*= (const root_of_unity &O)
            {
                std::rotate(p.rbegin(),p.rbegin()+O.position(s),p.rend());
                return reduce();
            }

            cyclic_ring operator*(const root_of_unity &O) const
            {
                cyclic_ring result=*this;
                return result*=O;
            }

            cyclic_ring& operator/= (const root_of_unity &O)
            {
                return (*this)*=(-O);
            }

            cyclic_ring operator/ (const root_of_unity &O) const
            {
                return (*this)*(-O);
            }

            cyclic_ring operator-(const cyclic_ring &O) const
            {
                cyclic_ring result;
                for(int i=0;i<s;i++)
                    result.p[i]=p[i]-O.p[i];
                return result;
            }

            cyclic_ring operator*(const cyclic_ring &O) const
            {
                cyclic_ring result;
                for(int i=0;i<s;i++) for(int j=0;j<s;j++)
                    result.p[(i+j)%s]+=p[i]*O.p[j];
                return result.reduce();
            }

            cyclic_ring& reduce()
            {
                auto deg=phi_s.degree();
                for(int i=s-1;i>=deg;i--) if(!is_zero(p[i]))
                {
                    auto c=p[i];
                    p[i]=0;
                    for(auto j:indexes)
                        p[i-deg+j]-=c*phi_s[j];
                }
                return *this;
            }

            cyclic_ring& operator+=(const cyclic_ring &O)
            {
                for(int i=0;i<s;i++)
                    p[i]+=O.p[i];
                return *this;
            }

            cyclic_ring& operator-=(const cyclic_ring &O)
            {
                for(int i=0;i<s;i++)
                    p[i]-=O.p[i];
                return *this;
            }

            cyclic_ring& operator*=(const cyclic_ring &O)
            {
                auto r=(*this)*O;
                p.swap(r.p);
                return *this;
            }


            cyclic_ring& operator+=(const root_of_unity &O)
            {
                auto m=O.position(s);
                p[m]+=1;
                return reduce();
            }

            cyclic_ring& operator-=(const root_of_unity &O)
            {
                auto m=O.position(s);
                p[m]-=1;
                return reduce();
            }

            cyclic_ring operator+(const cyclic_ring &O) const
            {
                cyclic_ring result=*this;
                return result+=O;
            }

            cyclic_ring operator-(const root_of_unity &O) const
            {
                cyclic_ring result=*this;
                return result-=O;
            }

            static void build_ring()
            {
                if(!phi_s.data().empty())
                    return;
                auto Phi= integer_cyclotomic_polynomial(s,F);
                phi_s.data().resize(Phi.data().size());
                for(int i=0;i<Phi.data().size()-1;i++) if(Phi[i]!=0)
                    indexes.push_back(i);
                Z_s.data().resize(s);
                Z_s.data().back()=1;
                Z_s[1]=-1;
                for(int i=0;i<Phi.data().size();i++)
                    phi_s.data()[i]=static_cast<R>(Phi.data()[i]);
                constexpr cp::integer M=998244353;
                using IK=cyclic<M>;
                polynomial<IK> Z,Phi_M;
                Phi_M.data().resize(Phi.data().size());
                for(int i=0;i<Phi.data().size();i++)
                    Phi_M.data()[i]=static_cast<IK>(Phi.data()[i]);
                Z.data().resize(s);
                Z.data().back()=1;
                Z[1]=-1;
                auto [A_M,B_M,G_M]= egcd(Z,Phi_M);
                auto g=G_M[0];
                auto K=g.inv()*s;
                for(auto &a:A_M)
                    a*=K;
                for(auto &b:B_M)
                    b*=K;
                for(auto &g:G_M)
                    g*=K;
                polynomial<cp::integer> A,B,G;
                if(G_M.degree()>0)
                    throw std::invalid_argument("Unable to build cyclic ring");
                A.data().resize(A_M.data().size());
                B.data().resize(B_M.data().size());
                G.data().resize(G_M.data().size());
                for(int i=0;i<A_M.data().size();i++)
                    A.data()[i]=static_cast<cp::integer>(A_M.data()[i]);
                for(int i=0;i<B_M.data().size();i++)
                    B.data()[i]=static_cast<cp::integer>(B_M.data()[i]);
                for(int i=0;i<G_M.data().size();i++)
                    G.data()[i]=static_cast<cp::integer>(G_M.data()[i]);
                for(auto &a:A) if(a>=M/2)
                        a-=M;
                for(auto &b:B) if(b>=M/2)
                        b-=M;
                for(auto &g:G) if(g>=M/2)
                        g-=M;
                H1.data().resize(A.data().size());
                H2.data().resize(B.data().size());
                for(int i=0;i<A.data().size();i++)
                    H1.data()[i]=static_cast<R>(A.data()[i]);
                for(int i=0;i<B.data().size();i++)
                    H2.data()[i]=static_cast<R>(B.data()[i]);
            }
            auto &operator[](int i)
            {
                return p[i];
            }
            const auto &operator[](int i) const
            {
                return p[i];
            }

        };

        template<typename R>
        struct cyclic_ring<R,2>
        {
            R p;
            cyclic_ring(const polynomial<R> &P):p{}
            {
                if(P.degree())
                    p=P[0];
            }
            cyclic_ring(const R &a=0):p{}
            {
                p=a;
            }

            template<std::integral A>
            cyclic_ring(const A &a):p{}
            {
                p=a;
            }

            cyclic_ring& operator*= (const root_of_unity &O)
            {
                if(O.position(2)%2)
                    p=-p;
                return reduce();
            }

            cyclic_ring operator*(const root_of_unity &O) const
            {
                cyclic_ring result=*this;
                return result*=O;
            }

            cyclic_ring& operator/= (const root_of_unity &O)
            {
                return (*this)*=O;
            }

            cyclic_ring operator/ (const root_of_unity &O) const
            {
                return (*this)*O;
            }

            cyclic_ring operator-(const cyclic_ring &O) const
            {
                return p-O.p;
            }

            cyclic_ring operator*(const cyclic_ring &O) const
            {
                return p*O.p;
            }

            cyclic_ring& reduce()
            {
                return *this;
            }

            cyclic_ring& operator+=(const cyclic_ring &O)
            {
                p+=O.p;
                return *this;
            }

            cyclic_ring& operator-=(const cyclic_ring &O)
            {
                p-=O.p;
                return *this;
            }

            cyclic_ring& operator*=(const cyclic_ring &O)
            {
                p*=O.p;
                return *this;
            }


            cyclic_ring& operator+=(const root_of_unity &O)
            {
                if(O.position(2)%2)
                    p-=1;
                else p+=1;
                return reduce();
            }

            cyclic_ring& operator-=(const root_of_unity &O)
            {
                if(O.position(2)%2)
                    p+=1;
                else p-=1;
                return reduce();
            }

            cyclic_ring operator+(const cyclic_ring &O) const
            {
                cyclic_ring result=*this;
                return result+=O;
            }

            cyclic_ring operator-(const root_of_unity &O) const
            {
                cyclic_ring result=*this;
                return result-=O;
            }

            inline static polynomial<R> H1,H2;
            static void build_ring()
            {
                H1={1};
                H2={1};
            }

            auto &operator[](int i)
            {
                return p;
            }
            const auto &operator[](int i) const
            {
                return p;
            }
        };

        template<typename R,int s>
        cyclic_ring<R,s> operator*(root_of_unity &a,const cyclic_ring<R,s> &b)
        {
            return b*a;
        }

        template<typename R,int s,bool inverse>
        struct cyclic_ring_extension;
        using namespace signals;

        template<typename T,int s>
        struct fixed_radix_fft : public abstract_fft<details::cyclic_ring_extension<T,s,false>>,
                                 public abstract_fft<details::cyclic_ring_extension<T,s,true>>
        {
            template<bool inverse>
            using R=cyclic_ring_extension<T,s,inverse>;
            using abstract_fft<R<false>>::transform;
            using abstract_fft<R<true>>::transform;

            template<bool id>
            void transform_rec(linalg::tensor_view<R<id>, 1> &v, bool inverse = false, FFTNormalization normalization = FFTNormalization::None) const {
                auto n = v.size();
                if (n == 1)
                    return;
                if (n % s != 0)
                    throw std::invalid_argument("size must be a multiple of s");
                std::size_t p = s, q = n / p;
                std::vector<linalg::tensor_subview<R<id>, 1>> V;
                for (unsigned i = 0; i < p; i++)
                    V.push_back(v.slice({i}, {n}, {p}));
                for (auto &v: V)
                    transform_rec(v, inverse, normalization);
                auto o=1;
                root_of_unity w(o,n),z(o*q,n);
                if (inverse) {
                    w = -w;
                    z = -z;
                }
                root_of_unity t(0,n);
                std::vector<R<id>> result(n);
                for (int i = 0; i < p; i++, t *= z) {
                    root_of_unity h1(0,n),h2(0,n);
                    for (int j = 0; j < p; j++, h1 *= t, h2 *= w) {
                        root_of_unity h3(0,n);
                        for (int k = 0; k < q; k++, h3 *= h2)
                            result[i * q + k] += h1 * h3 * V[j](k);
                    }
                }
                for (int i = 0; i < n; i++)
                    v(i) = result[i];
            }

            void transform(linalg::tensor_view<R<false>,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
            {
                transform_rec(v,inverse,normalization);
//                normalize(v,normalization);
            }

            void transform(linalg::tensor_view<R<true>,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
            {
                transform_rec(v,inverse,normalization);
//                normalize(v,normalization);
            }
        };

        template<typename T,int s> requires (std::has_single_bit<unsigned>(s))
        struct fixed_radix_fft<T,s> : public abstract_fft<details::cyclic_ring_extension<T,s,false>>,
                                 public abstract_fft<details::cyclic_ring_extension<T,s,true>>
        {
            template<bool inverse>
            using R = cyclic_ring_extension<T, s, inverse>;
            using abstract_fft<R<false>>::transform;
            using abstract_fft<R<true>>::transform;

            template<bool id>
            void transform_rec(linalg::tensor_view<R<id>, 1> &a, bool inverse = false,
                               FFTNormalization normalization = FFTNormalization::None) const {
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
                for (int len = 2; len <= n; len <<= 1)
                {
                    root_of_unity wlen(1,len);
                    if(inverse)
                        wlen=-wlen;
                    for (int i = 0; i < n; i += len)
                    {
                        root_of_unity w(0,n);
                        for (int j = 0; j < len / 2; j++)
                        {
                            R<id> u = a(i+j), v = w*a(i+j+len/2);
                            a(i+j) = u + v;
                            a(i+j+len/2) = u - v;
                            w *= wlen;
                        }
                    }
                }
            }
            void transform(linalg::tensor_view<R<false>,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
            {
                transform_rec(v,inverse,normalization);
//                normalize(v,normalization);
            }

            void transform(linalg::tensor_view<R<true>,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
            {
                transform_rec(v,inverse,normalization);
//                normalize(v,normalization);
            }
        };

        template<typename R,int s,bool inverse>
        void copy_to(const std::vector<cyclic_ring<R,s>> &X, std::vector<cyclic_ring_extension<R,s,inverse>> &Y,int l,int m)
        {
            Y.resize(l);
            for(int i=0;i<l;i++)
            {
                root_of_unity w(i,s*l);
                if(inverse)
                    w^=s-1;
                Y[i].data().resize(m);
                for(int j=0;j<m;j++)
                    Y[i][j]=X[i*m+j];
                Y[i]*=w;
            }
        }


        template<typename R,int r>
        void printf(std::ostream& H,const std::vector<cp::details::cyclic_ring<R,r>> &A)
        {
            for(auto &s:A)
            {
                std::cout << "[ ";
                if constexpr(r>2)
                    for(auto &t:s.p)
                        std::cout << t << ", ";
                else
                    std::cout << s.p << ", ";
                std::cout << "]\n";
            }

        }
        template<typename R,int r,bool inverse>
        void printf(std::ostream& H,const cp::details::cyclic_ring_extension<R,r,inverse>&s)
        {
            H << "[ ";
            for(auto a:s.data())
            {
                H << "( ";
                if constexpr (r>2)
                    for(auto b:a.p)
                        H << b << ", ";
                else
                    H << a.p << ", ";
                H << "); ";
            }
            H << "]\n";

        }


        template<typename R,int r,bool inverse>
        void printf(std::ostream& H,const std::vector<cp::details::cyclic_ring_extension<R,r,inverse>> &A)
        {
            for(auto &s:A)
            {
                H << "[ ";
                for(auto a:s.data())
                {
                    H << "( ";
                    if constexpr(r>2)
                        for(auto b:a.p)
                            H << b << ", ";
                    else
                        H << a.p << ", ";
                    H << "); ";
                }
                H << "]\n";
            }

        }

        template<typename R>
        std::vector<cyclic_ring<R,2>> cyclic_general_fast_multiplication_p2(const std::vector<cyclic_ring<R,2>> &P,const std::vector<cyclic_ring<R,2>> &Q,int n,
                                                                                        int threshold,int *mu = nullptr)
        {
            constexpr int h=1,K=1<<h;
            int k=std::bit_width<unsigned>(n)-1;
            using namespace details;
            using cyclic_ring=cyclic_ring<R,2>;
            if(n<=threshold)
            {
                cyclic_ring_extension<R,2,false> U,V,W;
                U.p=P;
                V.p=Q;
                W=U*V;
#ifndef NDEBUG
                std::cout << "--------Threshold Input-----------\n";
                printf(std::cout,U);
                std::cout << "------------------------\n";
                std::cout << "--------Threshold Output-----------\n";
                printf(std::cout,W);
                std::cout << "------------------------\n";
#endif
                return W.p;
            }
            int r=k/2;
            int l=1<<r,p=n/l,m=p<<1;
            using Ring=cyclic_ring_extension<R,2,false>;
            std::vector<Ring> A,B,C;
            A.resize(K*l);
            B.resize(K*l);
            for(int i=0;i<K*l;i++)
            {
                root_of_unity w(i,2*l);
                A[i].data().resize(m);
                B[i].data().resize(m);
                if(i<l)
                for(int j=0;j<p;j++)
                {
                    A[i][j]=P[i*p+j];
                    B[i][j]=Q[i*p+j];
                }
//                A[i]*=w;
//                B[i]*=w;
            }
            C.resize(K*l);
            static fixed_radix_fft<R,2> fft;
#ifndef NDEBUG
            std::cout << "--------Input-----------\n";
            printf(std::cout,A);
            std::cout << "------------------------\n";
#endif
            fft.transform(linalg::vector_view(A));
            fft.transform(linalg::vector_view(B));
#ifndef NDEBUG
            std::cout << "--------FFT-----------\n";
            printf(std::cout,A);
            std::cout << "------------------------\n";
#endif
            if(mu)
                (*mu)+=r+h;
            for(int i=0;i<K*l;i++)
                C[i].data()=cyclic_general_fast_multiplication_p2(A[i].data(),B[i].data(),m,threshold,i==0?mu:nullptr);
            for(int i=0;i<K*l;i++)
            {
//                root_of_unity w(i,2*l);
//                C[i]/=w;
            }
#ifndef NDEBUG
            std::cout << "--------Output-----------\n";
            printf(std::cout,C);
            std::cout << "------------------------\n";
#endif
            fft.transform(linalg::vector_view(C),true);
#ifndef NDEBUG
            std::cout << "--------Inverse FFT-----------\n";
            printf(std::cout,C);
            std::cout << "------------------------\n";
#endif
            std::vector<cyclic_ring> result(n);
            for(int i=0;i<K*l;i++) for(int j=0;j<m;j++)
                result[(i*p+j)%n]+=C[i][j]*root_of_unity((i*p+j)/n,2);
#ifndef NDEBUG
            std::cout << "--------Result-----------\n";
            printf(std::cout,result);
            std::cout << "------------------------\n";
#endif
            return result;
        }

        template<typename R>
        std::vector<cyclic_ring<R,2>> cyclic_general_fast_pow_p2(const std::vector<cyclic_ring<R,2>> &P, unsigned long long s,int n,
                                                                            int threshold,int *mu = nullptr)
        {
            int k=std::bit_width<unsigned>(n);
            using namespace details;
            using cyclic_ring=cyclic_ring<R,2>;
            cyclic_ring::build_ring();
            if(n<=threshold)
            {
                cyclic_ring_extension<R,2,false> U;
                U.p=P;
//                printf(std::cout,U);
                return pow(U,s).p;
            }
            int l=1<<(k/2),p=1<<((k+1)/2),m=p<<1;
            using Ring=cyclic_ring_extension<R,2,false>;
            std::vector<Ring> A,C;
            A.resize(l);
            for(int i=0;i<l;i++)
            {
                root_of_unity w(i,2*l);
                A[i].data().resize(m);
                for(int j=0;j<p;j++)
                    A[i][j]=P[i*p+j];
//                A[i]*=w;
//                B[i]*=w;
            }
            C.resize(l);
            static fixed_radix_fft<R,2> fft;
            fft.transform(linalg::vector_view(A));
            if(mu)
                (*mu)+=k/2+1;
            for(int i=0;i<l;i++)
                C[i].data()=cyclic_general_fast_pow_p2(A[i].data(),s,m,threshold,i==0?mu:nullptr);
            fft.transform(linalg::vector_view(C),true);
//            for(int i=0;i<l;i++)
//            {
//                root_of_unity w(i,2*l);
//                C[i]/=w;
//            }
            std::vector<cyclic_ring> result(n);
            for(int i=0;i<l;i++) for(int j=0;j<m;j++)
                    result[(i*p+j)%n]+=C[i][j]*root_of_unity(i*p+j >= n,2);
            return result;
        }


        template<bool id,int s,typename R>
        std::vector<cyclic_ring<R,s>> non_invertible_cyclic_general_fast_multiplication(const std::vector<cyclic_ring<R,s>> &P,const std::vector<cyclic_ring<R,s>> &Q,int n,int k,
                                                                                        int threshold,int *mu = nullptr)
        {
            using namespace details;
            using cyclic_ring=cyclic_ring<R,s>;
            cyclic_ring::build_ring();
            if(n<=threshold)
            {
                cyclic_ring_extension<R,s,id> P1,Q1;
                P1.p=P;
                Q1.p=Q;
//            printf(std::cout,P1);
//            std::cout << "-------------------\n";
//            printf(std::cout,Q1);
//            std::cout << "-------------------\n";
                return (P1*Q1).p;
            }
            int l=pow(s,k/2),m=n/l;
            using R1=cyclic_ring_extension<R,s,false>;
            using R2=cyclic_ring_extension<R,s,true>;
            std::pair<std::vector<R1>,std::vector<R2>> A,B,C;
            auto &[A1,A2]=A;
            auto &[B1,B2]=B;
            auto &[C1,C2]=C;
            copy_to(P,A1,l,m);
            copy_to(Q,B1,l,m);
            copy_to(P,A2,l,m);
            copy_to(Q,B2,l,m);
            C1.resize(l);
            C2.resize(l);
            fixed_radix_fft<R,s> fft;
//            printf(std::cout,A1);
//            std::cout << "-------------------\n";
//            printf(std::cout,A2);
//            std::cout << "-------------------\n";
//            printf(std::cout,B1);
//            std::cout << "-------------------\n";
//            printf(std::cout,B2);
//            std::cout << "-------------------\n";
            fft.transform(linalg::vector_view(A1));
            fft.transform(linalg::vector_view(A2));
            fft.transform(linalg::vector_view(B1));
            fft.transform(linalg::vector_view(B2));
//            std::cout << "---FFT-Transform---\n";
//            std::cout << "-------------------\n";
//
//            printf(std::cout,A1);
//            std::cout << "-------------------\n";
//            printf(std::cout,A2);
//            std::cout << "-------------------\n";
//            printf(std::cout,B1);
//            std::cout << "-------------------\n";
//            printf(std::cout,B2);
//            std::cout << "-------------------\n";
            if(mu)
                (*mu)+=k/2+1;
            for(int i=0;i<l;i++)
            {
                C1[i].data()=non_invertible_cyclic_general_fast_multiplication<false,s>(A1[i].data(),B1[i].data(),m,(k+1)/2,threshold,i==0?mu:nullptr);
                C2[i].data()=non_invertible_cyclic_general_fast_multiplication<true,s>(A2[i].data(),B2[i].data(),m,(k+1)/2,threshold,nullptr);
            }
//            printf(std::cout,C1);
//            std::cout << "-------------------\n";
//            printf(std::cout,C2);
//            std::cout << "-------------------\n";
//            std::cout << "-------------------\n";

            fft.transform(linalg::vector_view(C1),true);
            fft.transform(linalg::vector_view(C2),true);
            for(int i=0;i<l;i++)
            {
                root_of_unity w1(i,s*l),w2((s-1)*i,s*l);
                C1[i]/=w1;
                C2[i]/=w2;
            }
//            printf(std::cout,C1);
//            std::cout << "-------------------\n";
//            printf(std::cout,C2);
            cyclic_ring H(cyclic_ring::H1);
            std::vector<cyclic_ring> result(2*n);
            for(int i=0;i<l;i++)
            {
                auto dC=C2[i]-C1[i];
                auto dZ=H*dC;
                for(auto &c1:C1[i].data())
                    c1*=s;
                root_of_unity w(l,s*l);
                auto Z=C1[i]-dZ/w;
                for(int j=0;j<m;j++)
                {
                    result[i * m + j] += Z[j];
                    result[i * m + j + m] += dZ[j];
                }
            }
//            std::cout << "-------Result-------\n";
//            std::cout << "-------------------\n";
//            printf(std::cout,result);
//            std::cout << "-------------------\n";
            for(int i=n;i<2*n;i++)
            {
                root_of_unity w(id?s-1:1,s);
                result[i-n] += w*result[i];
            }
            result.resize(n);
//            std::cout << "-------Post--------\n";
//            std::cout << "-------------------\n";
//            printf(std::cout,result);
//            std::cout << "-------------------\n";
            return result;
        }

        template<typename R,int s,bool inverse>
        struct cyclic_ring_extension
        {
            std::vector<cyclic_ring<R,s>> p;
            cyclic_ring_extension(const std::vector<R> &P):p{}
            {
                for(int i=0;i<std::min<int>(s,P.size());i++)
                    p.emplace_back(P[i]);
            }

            cyclic_ring_extension(const R &a=0):p{}
            {
                p.emplace_back(a);
            }

            template<bool other>
            cyclic_ring_extension operator+(const cyclic_ring_extension<R,s,other> &O) const
            {
                cyclic_ring_extension result=*this;
                return result+=O;
            }

            template<bool other>
            cyclic_ring_extension operator-(const cyclic_ring_extension<R,s,other> &O) const
            {
                cyclic_ring_extension result=*this;
                return result-=O;
            }

            integer degree() const
            {
                return p.size();
            }

            auto &operator[](int i)
            {
                return p[i];
            }

            const auto &operator[](int i) const
            {
                return p[i];
            }

            auto& data()
            {
                return p;
            }
            const auto& data() const
            {
                return p;
            }

            template<bool other>
            cyclic_ring_extension& operator+=(const cyclic_ring_extension<R,s,other> &O)
            {
                if(O.degree()>degree())
                    p.resize(O.degree());
                for(int i=0;i<degree();i++)
                    p[i]+=O.p[i];
                return *this;
            }

            cyclic_ring_extension& operator+=(const root_of_unity &O)
            {
                auto m=O.position(s*degree());
                auto [q,r]=std::div(m,degree());
                p[r]+=root_of_unity(q,s);
                return *this;
            }

            cyclic_ring_extension& operator-=(const root_of_unity &O)
            {
                auto m=O.position(s*degree());
                auto [q,r]=std::div(m,degree());
                p[r]-=root_of_unity(q,s);
                return *this;
            }

            cyclic_ring_extension operator+(const root_of_unity &O) const
            {
                cyclic_ring_extension result=*this;
                return result+=O;
            }

            cyclic_ring_extension operator-(const root_of_unity &O) const
            {
                cyclic_ring_extension result=*this;
                return result-=O;
            }

            template<bool other>
            cyclic_ring_extension& operator-=(const cyclic_ring_extension<R,s,other> &O)
            {
                if(O.degree()>degree())
                    p.resize(O.degree());
                for(int i=0;i<degree();i++)
                    p[i]-=O.p[i];
                return *this;
            }

            cyclic_ring_extension& operator*=(const root_of_unity &O)
            {
                auto m=O.position(s*degree());
                auto [q,r]=std::div(m,degree());
                std::rotate(p.rbegin(),p.rbegin()+r,p.rend());
                for(int i=0;i<degree() && q+(i<r)>0;i++)
                {
                    auto w=root_of_unity(q + (i < r),s);
                    if constexpr (inverse)
                        p[i] /= w;
                    else
                        p[i] *= w;
                }
                return *this;
            }

            cyclic_ring_extension& operator/=(const root_of_unity &O)
            {
                auto m=O.position(s*degree());
                auto [q,r]=std::div(m,degree());
                std::rotate(p.begin(),p.begin()+r,p.end());
                for(int i=0;i<degree() && q+(i<r)>0;i++)
                {
                    auto w=root_of_unity(q + (i < r),s);
                    if constexpr (inverse)
                        p[p.size() - 1 - i] *= w;
                    else
                        p[p.size() - 1 - i] /= w;
                }
                return *this;
            }

            cyclic_ring_extension operator*(const root_of_unity &O) const
            {
                cyclic_ring_extension result=*this;
                return result*=O;
            }

            cyclic_ring_extension operator/(const root_of_unity &O) const
            {
                return (*this)*(-O);
            }

            cyclic_ring_extension operator*(const cyclic_ring_extension &b)
            {
                cyclic_ring_extension c;
                c.p.resize(p.size());
                auto m=p.size();
                for(int i=0;i<m;i++) for(int j=0;j<m;j++)
                {
                    auto [q,t]=std::div(i+j,m);
                    auto r=p[i]*b.p[j];
                    if constexpr (inverse)
                        r/=root_of_unity(q,s);
                    else
                        r*=root_of_unity(q,s);
                    c.p[t]+=r;
                }
                return c;
            }

            cyclic_ring_extension& operator*=(const cyclic_ring_extension &b)
            {
                return *this=*this*b;
            }
        };

        template<typename R,int s,bool inverse>
        cyclic_ring_extension<R,s,inverse> operator*(const root_of_unity &a,const cyclic_ring_extension<R,s,inverse> &b)
        {
            return b*a;
        }

        template<typename R,int s,bool inverse>
        cyclic_ring_extension<R,s,inverse> operator*(const cyclic_ring<R,s> &a,const cyclic_ring_extension<R,s,inverse> &B)
        {
            auto C=B;
            for(auto &c:C.data())
                c=a*c;
            return C;
        }
    }

    template<int s,typename R>
    std::pair<std::vector<R>,int> non_invertible_fast_general_polynomial_multiplication(const std::vector<R> &A,const std::vector<R> &B,int threshold=64)
    {
        using namespace details;
        std::vector<cyclic_ring<R,s>> U,V;
        for(auto &a:A)
            U.emplace_back(a);
        for(auto &b:B)
            V.emplace_back(b);
        auto n=A.size()+B.size()-1;
        int m=1,k=0;
        while(m<n)
        {
            m *= s;
            k++;
        }
        U.resize(m);
        V.resize(m);
        int mu=0;
        auto Z=non_invertible_cyclic_general_fast_multiplication<false,s>(U,V,m,k,std::max(threshold,s),&mu);
        std::vector<R> result(n);
        for(int i=0;i<n;i++)
            result[i]=Z[i][0];
        return std::make_pair(result,mu);
    }

    template<int s1=3,int s2=4,typename R>
    std::vector<R> fast_general_polynomial_multiplication(const std::vector<R> &A,const std::vector<R> &B,int threshold=64)
    {
        auto [U,mu1]=non_invertible_fast_general_polynomial_multiplication<s1>(A,B,threshold);
        auto [V,mu2]=non_invertible_fast_general_polynomial_multiplication<s2>(A,B,threshold);
        return {};
    }

    template<int s1=3,int s2=4,typename R>
    polynomial<R> fast_general_polynomial_multiplication(const polynomial<R> &A,const polynomial<R> &B,int threshold=64)
    {
        auto [U,mu1]=non_invertible_fast_general_polynomial_multiplication<s1>(A,B,threshold);
        auto [V,mu2]=non_invertible_fast_general_polynomial_multiplication<s2>(A,B,threshold);
        return {};
    }


    template<typename R>
    std::pair<std::vector<R>,int> fast_general_polynomial_multiplication_p2(const std::vector<R> &A,const std::vector<R> &B,int threshold=64)
    {
        using namespace details;
        std::vector<cyclic_ring<R,2>> U,V;
        for(auto &a:A)
            U.emplace_back(a);
        for(auto &b:B)
            V.emplace_back(b);
        auto n=A.size()+B.size()-1;
        int m=1,k=0;
        while(m<n)
        {
            m *= 2;
            k++;
        }
        U.resize(m);
        V.resize(m);
        int mu=0;
        auto Z=cyclic_general_fast_multiplication_p2(U,V,m,std::max(threshold,4),&mu);
        std::vector<R> result(n);
        for(int i=0;i<n;i++)
            result[i]=Z[i][0];
        return std::make_pair(result,mu);
    }

    template<typename R>
    std::pair<std::vector<R>,int> fast_general_polynomial_pow_p2(const std::vector<R> &A,unsigned long long s,int threshold=64)
    {
        using namespace details;
        std::vector<cyclic_ring<R,2>> U,V;
        for(auto &a:A)
            U.emplace_back(a);

        auto n=A.size();
        int m=1,k=0;
        while(m<n)
        {
            m *= 2;
            k++;
        }
        U.resize(m);
        V.resize(m);
        int mu=0;
        auto Z=cyclic_general_fast_pow_p2(U,s,m,std::max(threshold,8),&mu);
        std::vector<R> result(n);
        for(int i=0;i<n;i++)
            result[i]=Z[i][0];
        return std::make_pair(result,mu);
    }

}

#endif //CPLIBRARY_FAST_GENERAL_POLYNOMIAL_MULTIPLICATION_H

#include <iostream>
#include <chrono>
#include "nt/modular/static.h"
#include "nt/modular_functions.h"
#include "signals/ntt.h"
#include "fields/gf8.h"
#include "polynomial/fast_polynomial.h"
#include "nt/primality.h"
#include "signals/ntt.h"
#include "nt/modular/fixed.h"
#include "rings/quadratic/fixed.h"
#include "rings/quadratic/static.h"
#include "rings/quadratic/dynamic.h"




struct mixed_radix_fft : public cp::signals::abstract_fft<std::complex<double>>, protected cp::default_factoriser_t
{
    using R=std::complex<double>;
    using cp::signals::abstract_fft<R>::transform;
    std::shared_ptr<cp::abstract_factoriser> F;
    mixed_radix_fft(std::shared_ptr<cp::abstract_factoriser> _F=default_factoriser):F(_F){}
    void transform_rec(cp::linalg::tensor_view<R,1> &v, bool inverse=false, cp::signals::FFTNormalization normalization = cp::signals::FFTNormalization::None) const
    {
        auto n=v.size();
        if(n==1)
            return;
        std::uint32_t p = F->smallest_divisor(n);
        auto q=n/p;
        std::vector<cp::linalg::tensor_subview<R,1>> V;
        for(unsigned i=0;i<p;i++)
            V.push_back(v.slice({i},{n},{p}));
        for(auto &v:V)
            transform_rec(v,inverse,normalization);
        R w=std::polar(1.0,2*std::numbers::pi/n);
        if(inverse)
            w=std::conj(w);
        R z=std::pow(w,q);
        R t=1;
        std::vector<R> result(n);
        for(int i=0;i<p;i++,t*=z)
        {
            R h1=1,h2=1;
            for (int j = 0; j < p; j++,h1*=t,h2*=w)
            {
                R h3=1;
                for (int k = 0; k < q; k++,h3*=h2)
                    result[i*q+k] += h1 * h3 * V[j](k);
            }
        }
        for(int i=0;i<n;i++)
            v(i)=result[i];
    }

    void transform(cp::linalg::tensor_view<R,1> &v, bool inverse=false, cp::signals::FFTNormalization normalization = cp::signals::FFTNormalization::None) const override
    {
        transform_rec(v,inverse,normalization);
        cp::signals::normalize(v,normalization);
    }
};

/*
 * Fast Fourier Transform over Complex numbers
 * Suitable for smooth numbers
 * */
template<bool is_inverse=false>
struct fast_fourier
{
    using real=double;
    const real pi = acos(-1);
    int n;
    std::complex<real> w;
    using IC=std::complex<real>;
    inline static constexpr int sign=is_inverse?1:-1;
    inline static std::optional<std::reference_wrapper<cp::factoriser>> F_ref=std::optional<std::reference_wrapper<cp::factoriser>>();

public:
    inline static bool use_normalized=false;
    fast_fourier(int _n):n(_n),w(std::exp(IC(0,2*sign*pi/n)))
    {
    }
    virtual std::vector<IC> unnormalized(const std::vector<IC> &X) const
    {
        if(n==1)
            return X;
        auto &F=F_ref.value().get();
        auto p=F.smallest_divisor(n),q=n/p;
        fast_fourier<is_inverse> FFT(q);
        std::vector<std::vector<IC>> U(p,std::vector<IC>(q));
        for(int i=0;i<n;i++)
            U[i%p][i/p]=X[i];
        std::vector<std::vector<IC>> V(p);
        for(int i=0;i<p;i++)
            V[i]=FFT.unnormalized(U[i]);
        std::vector<IC> R(n);
        IC z=std::pow(w,q);
        IC t=1;
        for(int i=0;i<p;i++,t*=z)
        {
            IC h1=1,h2=1;
            for (int j = 0; j < p; j++,h1*=t,h2*=w)
            {
                IC h3=1;
                for (int k = 0; k < q; k++,h3*=h2)
                    R[i*q+k] += h1 * h3 * V[j][k];
            }
        }
        return R;
    }
    std::vector<IC> operator()(const std::vector<IC> &X) const
    {
        return use_normalized? normalized(X):unnormalized(X);
    }
    std::vector<IC> normalized(const std::vector<IC>&X) const
    {
        auto Y= unnormalized(X);
        for(auto &y:Y)
            y/=std::sqrt(n);
        return Y;
    }
    static void set_factoriser(cp::factoriser &F)
    {
        F_ref=F;
        fast_fourier<!is_inverse>::F_ref=F;
    }

    static cp::factoriser& get_factoriser()
    {
        return F_ref.value();
    }
};

void inplace_fft2(std::vector<std::complex<double>> & a, bool inverse=false) {
    int n = a.size();
    const double pi = acos(-1);
    using IC=std::complex<double>;
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;

        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * pi / len * (inverse ? -1 : 1);
        IC wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            IC w(1);
            for (int j = 0; j < len / 2; j++) {
                IC u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    if (inverse) {
        for (IC & x : a)
            x /= n;
    }
}

template<cp::integer m>
std::vector<cp::cyclic<m>> fast_multiplication_complex(const std::vector<cp::cyclic<m>> &a,const std::vector<cp::cyclic<m>> &b)
{
    using IC=std::complex<double>;
    using namespace cp;
    integer block=std::ceil(std::sqrt(m));
    integer n=std::bit_ceil(a.size()+b.size()-1);
    std::vector<IC> A(n),B(n);
    for(int i=0;i<a.size();i++)
    {
        auto [q,r]=std::div(static_cast<integer>(a[i]),block);
        A[i].real(r);
        A[i].imag(q);
    };
    for(int i=0;i<b.size();i++)
    {
        auto [q,r]=std::div(static_cast<integer>(b[i]),block);
        B[i].real(r);
        B[i].imag(q);
    };
    inplace_fft2(A);
    inplace_fft2(B);
    std::vector<IC> C(n),D(n);
    auto r=n;
    for(unsigned i=0;i<r;i++)
    {
        auto j=(r-i)&(r-1);
        C[i]=A[i]*A[i]-std::conj(A[j]*A[j]);
        C[i]*=IC(0,-0.25);
        D[i]=B[i]*B[j]-std::conj(B[j]*B[i]);
        D[i]*=IC(0,-0.25);
    }
    std::vector<IC> X(n),Y(n),Z(n);
    for(int i=0;i<n;i++)
    {
        X[i].real(C[i].real()*D[i].real());
        Y[i].imag(C[i].real()*D[i].imag()+C[i].imag()*D[i].real());
    }
    inplace_fft2(X,true);
    inplace_fft2(Y,true);
    std::vector<cyclic<m>> R(n);
    for(int i=0;i<n;i++)
    {
        R[i]+=std::llround(X[i].real());
        R[i]+=std::llround(X[i].imag())%m*block;
        R[i]+=std::llround(Y[i].real())%m*(block*block)%m;
    }
    return R;
}

constexpr cp::integer M=998244353;
int main()
{
    using R=std::complex<double>;
    using namespace cp;
    integer n;
    std::cin >> n;
    std::vector<R> B(n);
    for(int i=0;i<n;i++)
        B[i] = i+1;
    linalg::vector_view A(B);
    std::shared_ptr<cp::abstract_factoriser> F=std::make_shared<cp::factoriser>(1e6);
    signals::radix2_fft<R> FFT;
    auto t1=std::chrono::high_resolution_clock::now();
    FFT.transform(A);
    auto t2=std::chrono::high_resolution_clock::now();
}
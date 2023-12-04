//
// Created by ramizouari on 27/11/23.
//

#ifndef CPLIBRARY_FFT_H
#define CPLIBRARY_FFT_H
#include "linear_algebra/view.h"
#include <complex>
#include "algebra/bits.h"
#include <cstdint>
namespace cp::signals
{

    enum class FFTNormalization
    {
        None,
        Sqrt,
        Normalized,
        Auto
    };

    using cp::linalg::dynamic_extent;

    template<typename R>
    struct abstract_fft
    {
        virtual void transform(linalg::tensor_view<R,1> &&v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            transform(v,inverse,normalization);
        }
        virtual void transform(linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const = 0;
        virtual void transform(linalg::tensor_view<R,dynamic_extent> &&v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            transform(v,inverse,normalization);
        }
        virtual void transform(linalg::tensor_view<R,dynamic_extent> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            if(v.rank()>1)
                throw std::invalid_argument("rank must be one");
            auto x=to_static_view<1,R>(v);
            transform(x,inverse,normalization);
        }
        R& operator()(R &x, bool inverse=false, FFTNormalization normalization=FFTNormalization::None) const
        {
            transform(x,inverse,normalization);
            return x;
        }
        virtual ~abstract_fft()= default;
    };

    template<typename R>
    void normalize(linalg::tensor_view<std::complex<R>,1> &v,FFTNormalization normalized)
    {
        R r;
        switch (normalized)
        {
            case FFTNormalization::None:
                r=1;
                break;
            case FFTNormalization::Sqrt:
                r=std::sqrt(v.size());
                break;
            case FFTNormalization::Normalized:
                r=v.size();
                break;
            case FFTNormalization::Auto:
                throw std::invalid_argument("cannot normalize with auto");
        }
        if(normalized!=FFTNormalization::None) for (std::complex<R> & x : v)
            x /= r;
    }

    template<typename R>
    void normalize(linalg::tensor_view<std::complex<R>,1> &&v,FFTNormalization normalized)
    {
        normalize(v,normalized);
    }

    template<typename R>
    void inplace_fft2(cp::linalg::tensor_view<std::complex<R>,1> & a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        if(normalized==FFTNormalization::Auto)
            normalized=FFTNormalization::Sqrt;
        int n = a.size();
        for (int i = 1, j = 0; i < n; i++)
        {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1)
                j ^= bit;
            j ^= bit;
            if (i < j)
                swap(a(i), a(j));
        }
        for (int len = 2; len <= n; len <<= 1)
        {
            R ang = 2 * std::numbers::pi / len * (inverse ? -1 : 1);
            std::complex<R> wlen = std::polar<R>(1.,ang);
            for (int i = 0; i < n; i += len)
            {
                std::complex<R> w(1);
                for (int j = 0; j < len / 2; j++)
                {
                    std::complex<R> u = a(i+j), v = a(i+j+len/2) * w;
                    a(i+j) = u + v;
                    a(i+j+len/2) = u - v;
                    w *= wlen;
                }
            }
        }
        normalize(a,normalized);
    }

    template<typename R>
    void inplace_fft2(cp::linalg::tensor_view<std::complex<R>,1> && a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        inplace_fft2(a,inverse,normalized);
    }

    template<typename R>
    struct radix2_fft: public abstract_fft<R>
    {
        using abstract_fft<R>::transform;
        void transform(linalg::tensor_view<R,1> &v,bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
        {
            if(v.size()!= bit_ceil(v.size()))
                throw std::invalid_argument("size of vector must be a power of 2");
            inplace_fft2(v,inverse,normalization);
        }
    };

    template<typename R>
    struct mixed_radix_fft;


    template<std::floating_point Real>
    struct mixed_radix_fft<std::complex<Real>> : public cp::signals::abstract_fft<std::complex<Real>>, protected cp::default_factoriser_t
    {
        using R=std::complex<Real>;
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
            R z=std::polar(1.0,2*std::numbers::pi/p);
            if(inverse)
            {
                w = std::conj(w);
                z = std::conj(z);
            }
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
            normalize(v,normalization);
        }
    };

    template<typename R>
    R mod_operator(R x, R y)
    {
        return (x%y+y)%y;
    }

    //Bluestein's algorithm
    template<typename R>
    std::vector<std::complex<R>> general_fft(const cp::linalg::tensor_view<std::complex<R>,1> &a, bool inverse, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        using cp::linalg::vector_view;
        unsigned int n = a.size();
        if(bit_ceil(n)==n)
        {
            std::vector<std::complex<R>> b(n);
            for(int i=0;i<n;i++)
                b[i]=a(i);
            inplace_fft2(vector_view(b),inverse);
            return b;
        }
        //m>=2*n
        auto m=bit_ceil(2*n-1);
        std::vector<std::complex<R>> w(m);
        R ang = std::numbers::pi  / n * (inverse ? 1 : -1);
        w[1]=std::polar<R>(1.,ang);
        w[0]=1;
        for(int i=2;i<m;i++)
            w[i]=w[i-1]*w[1];
        std::vector<std::complex<R>> A(m),B(m),W(m);
        for(size_t i=0;i<n;i++)
        {
            auto r=mod_operator<std::int64_t>(i*i,2*n);
            W[i]=w[r];
            A[i]=a(i)*w[2*n-r];
            B[i]=w[r];
        }
        for(size_t i=1;i<n;i++)
            B[m-i]=B[i];
        inplace_fft2(vector_view(A),false,FFTNormalization::None);
        inplace_fft2(vector_view(B),false,FFTNormalization::None);
        for(size_t i=0;i<m;i++)
            A[i]*=B[i];
        inplace_fft2(vector_view(A),true,FFTNormalization::Normalized);
        for(size_t i=0;i<n;i++)
            A[i]*=std::conj(W[i]);
        A.resize(n);
        normalize(vector_view(A),normalized);
        return A;
    }

    template<typename R>
    struct bluestein_fft : public abstract_fft<R>
    {
        using abstract_fft<R>::transform;
        void transform(linalg::tensor_view<R,1> &v,bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
        {
            auto b=general_fft(v,inverse,normalization);
            for(int i=0;i<v.size();i++)
                v(i)=b[i];
        }
    };

    template<typename R>
    struct chirpz_transform
    {
        radix2_fft<R> fft;
        R z;
        template<typename ...Args>
        chirpz_transform(R _z,Args&&... args):z(_z),fft(std::forward<Args>(args)...)
        {
        }

        std::vector<R> transform(const cp::linalg::tensor_view<R,1> &a)
        {
            using linalg::vector_view;
            auto n=a.size();
            auto m=bit_ceil(2*n-1);
            std::vector<R> A(m),B(m),W(m);
            for(size_t i=0;i<n;i++)
            {
                W[i]=pow(z,i*i);
                A[i]=a(i)*W[i];
            }
            for(size_t i=0;i<n;i++)
                B[i]=R(1)/W[i];
            for(size_t i=1;i<n;i++)
                B[m-i]=B[i];
            fft.transform(vector_view(A),false,FFTNormalization::None);
            fft.transform(vector_view(B),false,FFTNormalization::None);
            for(size_t i=0;i<m;i++)
                A[i]*=B[i];
            fft.transform(vector_view(A),true,FFTNormalization::Normalized);
            for(size_t i=0;i<n;i++)
                A[i]/=W[i];
            A.resize(n);
            return A;
        }

        std::vector<R> transform(const cp::linalg::tensor_view<R,1> &&a)
        {
            return transform(a);
        }

    };

}

#endif //CPLIBRARY_FFT_H

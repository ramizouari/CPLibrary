//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_NTT_H
#define CPLIBRARY_NTT_H

#include <memory>
#include <set>
#include <optional>
#include <utility>
#include "fft.h"
#include "nt/modular_arithmetic.h"
#include "nt/modular_functions.h"
#include "nt/number_theory.h"

namespace cp::signals
{

    template<integer n>
    struct abstract_ntt : public abstract_fft<cyclic<n>>, protected default_factoriser_t
    {
        using R=cyclic<n>;
        using abstract_fft<R>::transform;
        explicit abstract_ntt(std::shared_ptr<abstract_factoriser> _F=default_factoriser):F(std::move(_F)){}
        mutable integer version=0;
        mutable std::array<std::unordered_map<integer,R>,2> cache;
        mutable std::optional<integer> phi;
        mutable R w1,w2;
        mutable std::shared_ptr<abstract_factoriser> F;

        void build() const
        {
            version=cyclic<n>::modulus();
            if(!F)
                F=default_factoriser;
            phi= carmichael_totient(cyclic<n>::modulus(),*F);
            w1= primitive_root_of_unity(cyclic<n>::modulus(),*F);
            w2=w1.pinv();
        }
        virtual R root_of_unity(integer size,integer m,bool inverse) const
        {
            if(version!=cyclic<n>::modulus())
                build();
            R w;
            if(cache[inverse].count(size))
                w=cache[inverse][size];
            else
            {
                w=inverse?w2:w1;
                auto [q,r]=std::div(*phi,size);
                if(r!=0)
                    throw std::invalid_argument("size must divide phi(m)");
                w=pow(w,q);
                cache[inverse][size]=w;
            }
            return w;
        }
        virtual ~abstract_ntt()= default;
    };

    template<>
    struct abstract_ntt<dynamic_modulus> : public abstract_fft<cyclic<dynamic_modulus>>, protected default_factoriser_t
    {
        using R=cyclic<dynamic_modulus>;
        using abstract_fft<R>::transform;
        explicit abstract_ntt(std::shared_ptr<abstract_factoriser> _F=default_factoriser):F(std::move(_F)){}
        mutable std::array<std::map<std::pair<integer,integer>,R>,2> cache;
        mutable std::map<integer,integer> phi;
        mutable std::map<integer,R> w1,w2;
        mutable std::set<integer> versions;
        mutable std::shared_ptr<abstract_factoriser> F;

        void build(integer m) const
        {
            versions.emplace(m);
            if(!F)
                F=default_factoriser;
            phi[m]= carmichael_totient(m,*F);
            w1[m]= primitive_root_of_unity(m,*F);
            w2[m]=w1[m].pinv();
        }
        virtual R root_of_unity(integer size,integer m,bool inverse) const
        {
            if(!versions.count(m))
                build(m);
            R w;
            if(cache[inverse].count({m,size}))
                w=cache[inverse][{m,size}];
            else
            {
                w=inverse?w2[m]:w1[m];
                auto [q,r]=std::div(phi[m],size);
                if(r!=0)
                    throw std::invalid_argument("size must divide phi(m)");
                w=pow(w,q);
                cache[inverse][{m,size}]=w;
            }
            return w;
        }
        virtual ~abstract_ntt()= default;
    };

    inline void cyclic_normalize(linalg::tensor_view<cyclic<dynamic_modulus>,1> &v,integer modulus,FFTNormalization normalized)
    {
        cyclic<dynamic_modulus> r(0,modulus);
        switch (normalized)
        {
            case FFTNormalization::None:
                r=1;
                break;
            case FFTNormalization::Sqrt:
                //r=cp::sqrt(cyclic<dynamic_modulus>(v.size(),modulus));
                throw std::invalid_argument("not implemented");
                break;
            case FFTNormalization::Normalized:
                r=v.size();
                break;
        }
        r=r.pinv();
        for(auto &x:v)
            x*=r;
    }

    template<integer m>
    void cyclic_normalize(linalg::tensor_view<cyclic<m>,1> &v,FFTNormalization normalized)
    {
        if constexpr (m==dynamic_modulus)
        {
            return cyclic_normalize(v,v(0).modulus(),normalized);
        }
        else
        {
            cyclic<m> r;
            switch (normalized)
            {
                case FFTNormalization::None:
                    r=1;
                    break;
                case FFTNormalization::Sqrt:
                    //r=cp::sqrt(cyclic<m>(v.size()));
                    throw std::invalid_argument("not implemented");
                    break;
                case FFTNormalization::Normalized:
                    r=v.size();
                    break;
            }
            r=r.pinv();
            for(auto &x:v)
                x*=r;
        }
    }

    template<integer m>
    void cyclic_normalize(linalg::tensor_view<cyclic<m>,1> &&v,FFTNormalization normalized)
    {
        cyclic_normalize(v,normalized);
    }



    template<integer m>
    void inplace_ntt2(cp::linalg::tensor_view<cyclic<m>,1> & a, cyclic<m> theta, FFTNormalization normalized = FFTNormalization::Sqrt)
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
        std::vector<cyclic<m>> W;
        W.push_back(theta);
        for (int len = 2; len <= n; len <<= 1)
            W.push_back(W.back()*W.back());
        W.pop_back();
        for (int len = 2,r=W.size()-1; len <= n; len <<= 1,r--)
        {
            auto wlen = W[r];
            for (int i = 0; i < n; i += len)
            {
                cyclic<m> w(1);
                for (int j = 0; j < len / 2; j++)
                {
                    cyclic<m> u = a(i+j), v = a(i+j+len/2) * w;
                    a(i+j) = u + v;
                    a(i+j+len/2) = u - v;
                    w *= wlen;
                }
            }
        }
        cyclic_normalize(a,normalized);
    }

    void inplace_ntt2(cp::linalg::tensor_view<cyclic<dynamic_modulus>,1> & a, cyclic<dynamic_modulus> theta, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        using Cyclic=cyclic<dynamic_modulus>;
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
        std::vector<Cyclic> W;
        W.push_back(theta);
        for (int len = 2; len <= n; len <<= 1)
            W.push_back(W.back()*W.back());
        W.pop_back();
        for (int len = 2,r=W.size()-1; len <= n; len <<= 1,r--)
        {
            auto wlen = W[r];
            for (int i = 0; i < n; i += len)
            {
                Cyclic w(1,theta.modulus());
                for (int j = 0; j < len / 2; j++)
                {
                    Cyclic u = a(i+j), v = a(i+j+len/2) * w;
                    a(i+j) = u + v;
                    a(i+j+len/2) = u - v;
                    w *= wlen;
                }
            }
        }
        cyclic_normalize(a,normalized);
    }

    template<typename Cyclic>
    void inplace_ntt2(cp::linalg::tensor_view<Cyclic,1> && a, Cyclic theta, FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        inplace_ntt2(a,theta,normalized);
    }

    template<integer n>
    struct radix2_fft<cyclic<n>> : public abstract_ntt<n>, protected default_factoriser_t
    {
        using R=cyclic<n>;
        using abstract_fft<cyclic<n>>::transform;
        using abstract_ntt<n>::root_of_unity;
        using abstract_ntt<n>::build;
        radix2_fft(std::shared_ptr<abstract_factoriser> _F=default_factoriser):abstract_ntt<n>(_F)
        {
        }
        void transform(linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
        {
            auto w= root_of_unity(v.size(),cyclic<n>::modulus(),inverse);
            if(normalization==FFTNormalization::Auto)
                normalization=inverse?FFTNormalization::Normalized:FFTNormalization::None;
            inplace_ntt2(v, w, normalization);
        }
    };

    template<>
    struct radix2_fft<cyclic<dynamic_modulus>> : public abstract_ntt<dynamic_modulus>
    {
        using R=cyclic<dynamic_modulus>;
        using abstract_fft::transform;
        using abstract_ntt::root_of_unity;
        using abstract_ntt::build;
        using default_factoriser_t::default_factoriser;
        radix2_fft(std::shared_ptr<abstract_factoriser> _F=default_factoriser):abstract_ntt(_F)
        {
        }
        void transform(linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
        {
            if(!std::has_single_bit(v.size()))
                throw std::invalid_argument("size must be a power of 2");
            auto w= root_of_unity(v.size(),v(0).modulus(),inverse);
            if(normalization==FFTNormalization::Auto)
                normalization=inverse?FFTNormalization::Normalized:FFTNormalization::None;
            inplace_ntt2(v, w, normalization);
        }
    };


    //Bluestein's algorithm
    template<integer q>
    std::vector<cyclic<q>> general_ntt(const abstract_ntt<q> &ntt,const cp::linalg::tensor_view<cyclic<q>,1> &a, bool inverse,
                                       FFTNormalization normalized = FFTNormalization::Sqrt)
    {
        using cp::linalg::vector_view;
        unsigned int n = a.size();
        using R=cyclic<q>;
        if(bit_ceil(n)==n)
        {
            std::vector<R> b(n);
            for(int i=0;i<n;i++)
                b[i]=a(i);
            ntt.transform(vector_view(b),inverse);
            cyclic_normalize(vector_view(b),normalized);
            return b;
        }
        //m>=2*n
        auto m=bit_ceil(2*n-1);
        std::vector<cyclic<q>> w(m);
        w[1]=ntt.root_of_unity(2*n,q,inverse);
        w[0]=1;
        for(int i=2;i<m;i++)
            w[i]=w[i-1]*w[1];
        std::vector<cyclic<q>> A(m),B(m),W1(m),W2(m);
        for(size_t i=0;i<n;i++)
        {
            auto r=mod_operator<std::int64_t>(i*i,2*n);
            W1[i]=w[r];
            W2[i]=w[2*n-r];
            A[i]=a(i)*w[2*n-r];
            B[i]=w[r];
        }
        for(size_t i=1;i<n;i++)
            B[m-i]=B[i];
        ntt.transform(vector_view(A),false,FFTNormalization::None);
        ntt.transform(vector_view(B),false,FFTNormalization::None);
        for(size_t i=0;i<m;i++)
            A[i]*=B[i];
        ntt.transform(vector_view(A),true,FFTNormalization::Normalized);
        for(size_t i=0;i<n;i++)
            A[i]*=W2[i];
        A.resize(n);
        cyclic_normalize(vector_view(A),normalized);
        return A;
    }

    template<integer n>
    struct bluestein_fft<cyclic<n>> : public abstract_fft<cyclic<n>>, private default_factoriser_t
    {
        using abstract_fft<cyclic<n>>::transform;
        using default_factoriser_t::default_factoriser;
        using R=cyclic<n>;
        std::unique_ptr<radix2_fft<R>> fft;
        bluestein_fft(std::shared_ptr<abstract_factoriser> _F=default_factoriser):fft(std::make_unique<radix2_fft<R>>(_F)){}
        void transform(linalg::tensor_view<R,1> &v,bool inverse, FFTNormalization normalization = FFTNormalization::None) const override
        {
            auto b=general_ntt(*fft,v,inverse,normalization);
            for(int i=0;i<v.size();i++)
                v(i)=b[i];
        }
    };

    template<integer m>
    struct mixed_radix_fft<cyclic<m>> : public abstract_ntt<m>, private default_factoriser_t
    {
        using R=cyclic<m>;
        using cp::signals::abstract_fft<R>::transform;
        using default_factoriser_t::default_factoriser;
        mixed_radix_fft(std::shared_ptr<cp::abstract_factoriser> _F=default_factoriser):abstract_ntt<m>(_F){}
        void transform_rec(cp::linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = cp::signals::FFTNormalization::None) const
        {
            auto n=v.size();
            if(n==1)
                return;
            std::uint32_t p = this->F->smallest_divisor(n);
            auto q=n/p;
            std::vector<cp::linalg::tensor_subview<R,1>> V;
            for(unsigned i=0;i<p;i++)
                V.push_back(v.slice({i},{n},{p}));
            for(auto &v:V)
                transform_rec(v,inverse,normalization);
            R w=this->root_of_unity(n,m,inverse);
            R z=this->root_of_unity(p,m,inverse);
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

        void transform(cp::linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = cp::signals::FFTNormalization::None) const override
        {
            transform_rec(v,inverse,normalization);
            cyclic_normalize(v,normalization);
        }
    };

}

#endif //CPLIBRARY_NTT_H

//
// Created by ramizouari on 27/11/23.
//

#ifndef CPLIBRARY_HADAMARD_H
#define CPLIBRARY_HADAMARD_H
#include "fft.h"
#include "multi_fft.h"
#include "algebra/bits.h"
namespace cp::signals
{

    template<typename R>
    void normalize_hadamard(linalg::tensor_view<R,1> &v,FFTNormalization normalized)
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
        }
        if(normalized!=FFTNormalization::None) for (R & x : v)
            x /= r;
    }
    template<typename R>
    struct binary_fft : public abstract_fft<R>
    {
        void transform(linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const override
        {
            if(v.size()>2)
                throw std::invalid_argument("size must be one or two");
            else if(v.size()==2)
            {
                R a=v(0),b=v(1);
                v(0)=a+b;
                v(1)=a-b;
            }
            normalize_hadamard(v,normalization);
        }
    };

    template<typename R>
    struct binary_view : public linalg::tensor_view<R,linalg::dynamic_extent>
    {
        linalg::tensor_view<R,1> &src;
        std::vector<std::size_t> m_shape;
        binary_view(linalg::tensor_view<R,1> &_src):src(_src)
        {
            if(src.size() != cp::bit_ceil(src.size()))
                throw std::invalid_argument("size must be a power of two");
            auto n=src.size();
            auto m=0;
            while(n>1)
            {
                n/=2;
                m++;
            }
            m_shape.resize(m,2);
        }
        R& at(std::vector<std::size_t> indexes)
        {
            std::size_t k=0;
            for(auto i:indexes)
                k=2*k+i;
            return src.at(k);
        }

        const R& at(std::vector<std::size_t> indexes) const
        {
            std::size_t k=0;
            for(auto i:indexes)
                k=2*k+i;
            return src.at(k);
        }

        std::size_t size() const
        {
            return src.size();
        }

        std::vector<std::size_t> shape() const
        {
            return m_shape;
        }
    };

    template<typename R>
    struct fast_hadamard
    {
        multi_fft<R> F;
    public:
        fast_hadamard():F(std::make_shared<binary_fft<R>>()){}
        void transform(linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            F.transform(binary_view<R>(v),inverse,normalization);
        }

        void transform(linalg::tensor_view<R,1> &&v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            F.transform(binary_view<R>(v),inverse,normalization);
        }
    };

    template<typename R>
    struct faster_hadamard
    {
        void transform(linalg::tensor_view<R,1> &v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            std::size_t h = 1;
            while (h < v.size())
            {
                for (std::size_t i=0;i< v.size(); i+=h * 2) for (std::size_t j=i; j< i + h;j++)
                {
                    auto x = v(j);
                    auto y = v(j + h);
                    v(j) = x + y;
                    v(j + h) = x - y;
                }
                h *= 2;
            }
            normalize_hadamard(v,normalization);

        }

        void transform(linalg::tensor_view<R,1> &&v, bool inverse=false, FFTNormalization normalization = FFTNormalization::None) const
        {
            transform(v,inverse,normalization);
        }
    };
}
#endif //CPLIBRARY_HADAMARD_H

//
// Created by ramizouari on 27/11/23.
//

#ifndef CPLIBRARY_MULTI_FFT_H
#define CPLIBRARY_MULTI_FFT_H

#include <memory>
#include <algorithm>
#include <utility>
#include "fft.h"
#include "linear_algebra/tensor.h"
#include "linear_algebra/utils.h"

namespace cp::signals
{
    template<typename R, size_t Rank>
    struct tensor_projection_view : linalg::tensor_view<R,1>
    {
        size_t k;
        linalg::tensor_view<R,Rank>& src;
        std::array<std::size_t,Rank> fixed;
        tensor_projection_view(linalg::tensor_view<R,Rank>& src,size_t k,
                             std::array<std::size_t,Rank> fixed):src(src),k(k),fixed(fixed){}
        R& at(std::array<std::size_t,1> indexes)
        {
            return at(indexes[0]);
        }

        R& at(size_t r)
        {
            auto H=fixed;
            H[k]=r;
            return src.at(H);
        }

        const R& at(size_t r) const
        {
            auto H=fixed;
            H[k]=r;
            return src.at(H);
        }

        const R& at(std::array<std::size_t,1> indexes) const
        {
            return at(indexes[0]);
        }

        std::array<std::size_t,1> shape() const
        {
            return {src.shape()[k]};
        }
    };

    template<typename R>
    struct tensor_projection_view<R,dynamic_extent> : public cp::linalg::tensor_view<R,1>
    {
        size_t k;
        linalg::tensor_view<R,dynamic_extent>& src;
        std::vector<std::size_t> fixed;
        tensor_projection_view(linalg::tensor_view<R,dynamic_extent>& src,size_t k,
                               std::vector<std::size_t> fixed):src(src),k(k),fixed(std::move(fixed)){}
        R& at(std::array<std::size_t,1> indexes)
        {
            return at(indexes[0]);
        }

        R& at(size_t r)
        {
            auto H=fixed;
            H[k]=r;
            return src.at(H);
        }

        const R& at(size_t r) const
        {
            auto H=fixed;
            H[k]=r;
            return src.at(H);
        }

        const R& at(std::array<std::size_t,1> indexes) const
        {
            return at(indexes[0]);
        }

        std::array<std::size_t,1> shape() const
        {
            return {src.shape()[k]};
        }
    };

    template<typename R>
    struct multi_fft:abstract_fft<R>
    {
        std::shared_ptr<abstract_fft<R>> F;
        multi_fft(std::shared_ptr<abstract_fft<R>> F):F(F)
        {
        }

        template<size_t Rank>
        void transform(linalg::tensor_view<R,Rank> &A, bool inverse=false, FFTNormalization normalized = FFTNormalization::Sqrt) const
        {
            auto shape=A.shape();
            for(int k=0;k<A.rank();k++)
            {
                std::array<std::size_t,Rank> I{},S=shape;
                S[k]=1;
                do {
                    tensor_projection_view<R,Rank> B(A,k,I);
                    linalg::flat_tensor<R,1> P({shape[k]}); // It is faster to use a new tensor due to cache
                    for(int i=0;i < shape[k];i++)
                        P(i)=B(i);
                    F->transform(P,inverse,normalized);
                    for(int i=0;i < shape[k];i++)
                        B(i)=P(i);
                    increment(I,S);
                } while(std::any_of(I.begin(),I.end(),[&](int i){return i!=0;}));
            }
        }

        template<size_t Rank>
        void transform(linalg::tensor_view<R,Rank> &&A, bool inverse=false, FFTNormalization normalized = FFTNormalization::Sqrt) const
        {
            transform(A,inverse,normalized);
        }

        void transform(linalg::tensor_view<R,dynamic_extent> &A, bool inverse=false, FFTNormalization normalized = FFTNormalization::Sqrt) const
        {
            auto shape=A.shape();
            for(int k=0;k<A.rank();k++)
            {
                std::vector<std::size_t> I(A.rank()),S=shape;
                S[k]=1;
                do {
                    tensor_projection_view<R,dynamic_extent> B(A,k,I);
                    linalg::flat_tensor<R,1> P({shape[k]}); // It is faster to use a new tensor due to cache
                    for(int i=0;i < shape[k];i++)
                        P(i)=B(i);
                    F->transform(P,inverse,normalized);
                    for(int i=0;i < shape[k];i++)
                        B(i)=P(i);
                    increment(I,S);
                }while(std::any_of(I.begin(),I.end(),[&](int i){return i!=0;}));
            }
        }

        void transform(linalg::tensor_view<R,1> &v, bool inverse, FFTNormalization normalization = FFTNormalization::Sqrt) const override
        {
            F->transform(v,inverse,normalization);
        }
    };
}

#endif //CPLIBRARY_MULTI_FFT_H

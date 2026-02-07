//
// Created by ramizouari on 14/12/24.
//

#ifndef CP_LIBRARY_MATRIX_OPERATIONS_H
#define CP_LIBRARY_MATRIX_OPERATIONS_H
#include "../../tensors/view.h"
#include "algebra/structures.h"
#include "../../tensors/tensor.h"

namespace cp::linalg {
    constexpr std::size_t matrix_output_rank(std::size_t rnk1,std::size_t rnk2) {
        return std::max(rnk1,rnk2)==dynamic_extent?dynamic_extent:rnk1+rnk2-2;
    }

    constexpr std::size_t einstein_output_rank(std::size_t rnk1,std::size_t rnk2,std::size_t dims) {
        return std::max(rnk1,rnk2)==dynamic_extent?dynamic_extent:rnk1+rnk2-2*dims;
    }


    template<additive_monoid M,std::size_t Rnk1,std::size_t Rnk2>
    void matrix_multiply(tensor_view<M,matrix_output_rank(Rnk1,Rnk2)> &c,const tensor_view<M,Rnk1> &a,const tensor_view<M,Rnk2> &b) {
        using OutType=tensor_view<M,matrix_output_rank(Rnk1,Rnk2)>;
        using It= typename OutType::iterator;
        auto I = a.begin().indexes;
        auto J = b.begin().indexes;
        auto r = b.shape()[0];
        for (It it = c.begin(); it != c.end(); ++it)
        {
            auto K = it.indexes;
            std::copy(K.begin(),K.end()-1,I.begin());
            std::copy(K.begin()+1,K.end(),J.begin());
            auto & v = *it;
            for (int k=0;k<r;k++) {
                I.back()=k;
                J.front()=k;
                v += a.at(I) * b.at(J);
            }
        }
    }

    template<additive_monoid M,std::size_t Rnk1,std::size_t Rnk2>
    flat_tensor<M,matrix_output_rank(Rnk1,Rnk2)> matrix_multiply(const tensor_view<M,Rnk1> &a,const tensor_view<M,Rnk2> &b)
    {
        constexpr auto Rnk3=matrix_output_rank(Rnk1,Rnk2);
        typename flat_tensor<M,Rnk3>::index_array outShape;
        auto in1Shape = a.shape();
        auto in2Shape = b.shape();
        if constexpr (Rnk3 == dynamic_extent) outShape.resize(in1Shape.size()+in2Shape.size()-2);
        std::copy(in1Shape.begin(),in1Shape.end()-1,outShape.begin());
        std::copy(in2Shape.begin()+1,in2Shape.end(),outShape.begin()+a.rank()-1);
        flat_tensor<M,Rnk3> c(outShape);
        matrix_multiply(c,a,b);
        return c;
    }


    template<std::size_t dims,additive_monoid M,std::size_t Rnk1,std::size_t Rnk2>
        void einstein_multiply(tensor_view<M,einstein_output_rank(Rnk1,Rnk2,dims)> &c,const tensor_view<M,Rnk1> &a,const tensor_view<M,Rnk2> &b)
    {
        auto I = a.begin().indexes;
        auto J = b.begin().indexes;
        std::array<std::size_t,dims> L;
        auto bShape = b.shape();
        std::copy(bShape.begin(),bShape.begin()+dims,L.begin());
        auto S = c.shape(),K = S;
        std::fill(K.begin(),K.end(),0);
        auto R = L;
        do {
            std::copy(K.begin(),K.end()-dims,I.begin());
            std::copy(K.begin()+dims,K.end(),J.begin()+dims);
            auto &v = c.at(K);
            do {
                std::copy(R.begin(),R.end(),I.rbegin());
                std::copy(R.begin(),R.end(),J.begin());
                v+=a.at(I) * b.at(J);
            }while (!rincrement(R,L));
        } while (!rincrement(K,S));
    }

    template<std::size_t dims,additive_monoid M,std::size_t Rnk1,std::size_t Rnk2>
        void einstein_multiply_transpose(tensor_view<M,einstein_output_rank(Rnk1,Rnk2,dims)> &c,const tensor_view<M,Rnk1> &a,const tensor_view<M,Rnk2> &b)
    {
        auto I = a.begin().indexes;
        auto J = b.begin().indexes;
        std::array<std::size_t,dims> L;
        auto bShape = b.shape();
        std::copy(bShape.end()-dims,bShape.end(),L.begin());
        auto S = c.shape(),K = S;
        std::fill(K.begin(),K.end(),0);
        auto R = L;
        std::fill(R.begin(),R.end(),0);
        do {
            std::copy(K.begin(),K.begin()+a.rank()-dims,I.begin());
            std::copy(K.begin()+a.rank()-dims,K.end(),J.begin());
            auto &v = c.at(K);
            do {
                std::copy(R.begin(),R.end(),I.end()-dims);
                std::copy(R.begin(),R.end(),J.end()-dims);
                v+=a.at(I) * b.at(J);
            }while (!rincrement(R,L));
        } while (!rincrement(K,S));
    }


    template<std::size_t dims,additive_monoid M,std::size_t Rnk1,std::size_t Rnk2>
    flat_tensor<M,einstein_output_rank(Rnk1,Rnk2,dims)> einstein_multiply(const tensor_view<M,Rnk1> &a,const tensor_view<M,Rnk2> &b)
    {
        constexpr auto Rnk3=einstein_output_rank(Rnk1,Rnk2,dims);
        typename flat_tensor<M,Rnk3>::index_array outShape;
        auto in1Shape = a.shape();
        auto in2Shape = b.shape();
        if constexpr (Rnk3 == dynamic_extent) outShape.resize(in1Shape.size()+in2Shape.size()-2*dims);
        std::copy(in1Shape.begin(),in1Shape.end()-dims,outShape.begin());
        std::copy(in2Shape.begin()+dims,in2Shape.end(),outShape.begin()+a.rank()-dims);
        flat_tensor<M,Rnk3> c(outShape);
        einstein_multiply<dims>(c,a,b);
        return c;
    }

    template<std::size_t dims,additive_monoid M,std::size_t Rnk1,std::size_t Rnk2>
    flat_tensor<M,einstein_output_rank(Rnk1,Rnk2,dims)> einstein_multiply_transpose(const tensor_view<M,Rnk1> &a,const tensor_view<M,Rnk2> &b)
    {
        constexpr auto Rnk3=einstein_output_rank(Rnk1,Rnk2,dims);
        typename flat_tensor<M,Rnk3>::index_array outShape;
        auto in1Shape = a.shape();
        auto in2Shape = b.shape();
        if constexpr (Rnk3 == dynamic_extent) outShape.resize(in1Shape.size()+in2Shape.size()-2*dims);
        std::copy(in1Shape.begin(),in1Shape.end()-dims,outShape.begin());
        std::copy(in2Shape.begin(),in2Shape.end()-dims,outShape.begin()+a.rank()-dims);
        flat_tensor<M,Rnk3> c(outShape);
        einstein_multiply_transpose<dims>(c,a,b);
        return c;
    }

}

#endif //CP_LIBRARY_MATRIX_OPERATIONS_H

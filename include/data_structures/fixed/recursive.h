//
// Created by ramizouari on 09/01/2026.
//

#ifndef CPLIBRARY_RECURSIVE_H
#define CPLIBRARY_RECURSIVE_H
#include <concepts>
#include <vector>

#include "segment_tree.h"
#include "algebra/binary_operation.h"
#include "functional/zip.h"

namespace cp::data_structures::fixed
{
    template<template <typename> typename RQ, typename O>
    concept range_query = requires {
        typename RQ<O>::key_type;
        typename RQ<O>::value_type;
        typename RQ<O>::binary_operation;
    } && requires(typename RQ<O>::key_type a, typename RQ<O>::key_type b) {
        true;//{ query(a, b) } -> std::convertible_to<typename RQ<O>::value_type>;
    };

    template<template <typename> typename RQ, typename O>
    concept supports_view = requires(RQ<O> rq) {
        typename O::value_type;
        {rq.span()} -> std::convertible_to<std::span<typename O::value_type>>;
    };

    template<template <typename> typename RQ, typename O> requires range_query<RQ,O>
    auto get_view_or_data(const RQ<O> &rq)
    {
        if constexpr (supports_view<RQ,O>)
            return rq.span();
        else return rq.data();
    }

    template<template <typename> typename RQ, typename O>
    concept range_query_update = range_query<RQ,O> && requires(typename RQ<O>::key_type a,
        typename RQ<O>::key_type b,
        typename RQ<O>::value_type c) {
        update(a, b, c);
    };

    template<template <typename> typename RQ,typename T, std::size_t Rank> requires range_query<RQ,T>
    struct recursive_range_query;

    template<template <typename> typename RQ, typename O, std::size_t Rank> requires range_query<RQ,O>
    struct recursive_range_query_operator : binary_operation<std::optional<recursive_range_query<RQ,O, Rank>>>
    {
        using binary_operation = recursive_range_query_operator<RQ,O, Rank - 1>;
        using T = std::optional<recursive_range_query<RQ,O, Rank>>;
        inline static recursive_range_query_operator<RQ,O, Rank - 1> op{};

        inline static std::optional<recursive_range_query<RQ,O, Rank>> neutral{};

        T reduce(const T& A, const T& B) const override
        {
            if (!A.has_value())
                return B;
            if (!B.has_value())
                return A;
            std::vector<std::optional<recursive_range_query<RQ,O, Rank - 1>>> result;
            auto A_ = get_view_or_data(*A);
            auto B_ = get_view_or_data(*B);
            for (auto [a,b]: zip(A_,B_))
                result.push_back(op(a,b));
            return result;
        }


    };

    template<template <typename> typename RQ, typename O> requires range_query<RQ,O>
    struct recursive_range_query_operator<RQ,O,0> : binary_operation<std::optional<recursive_range_query<RQ,O,0>>>
    {
        using T=std::optional<recursive_range_query<RQ,O,0>>;
        inline static O op{};
        inline static T neutral{};

        T reduce(const T& A, const T& B) const override
        {
            if (!A.has_value())
                return B;
            if (!B.has_value())
                return A;
            return op(*A,*B);
        }
    };

    template<typename Src, typename Target=Src>
    auto make_vector_optional(const std::vector<Src> &v)
    {
        std::vector<std::optional<Target>> res;
        res.reserve(v.size());
        for (auto &x:v)
            res.push_back(std::optional<Target>(std::in_place, x));
        return res;
    }

    template<template <typename> typename RQ,typename O, std::size_t Rank> requires range_query<RQ,O>
    struct recursive_range_query : protected RQ<recursive_range_query_operator<RQ,O, Rank - 1>>
    {
        using base_type = RQ<recursive_range_query_operator<RQ, O, Rank - 1>>;
        using value_type = RQ<O>::value_type;
        using key_type = std::array<int,Rank>;
        using key_view_type = std::span<const int,Rank>;
        using child_type = recursive_range_query<RQ, O, Rank - 1>;
        using tensor_type = std::vector<typename child_type::tensor_type>;

        friend class recursive_range_query_operator<RQ,O, Rank>;
        public:
        using base_type::base_type;

        recursive_range_query() = default;

        recursive_range_query(const tensor_type &elements) :
            base_type(std::vector(std::from_range,
                make_vector_optional<typename child_type::tensor_type,child_type>(elements)))
        {
        }

        template<typename RangeU, typename RangeV>
        value_type query(RangeU U, RangeV V)
        {
            key_view_type A(U),B(V);
            return query(A,B);
        }

        value_type query(key_view_type U,key_view_type V) = delete;


    };

    template<typename O, std::size_t Rank>
    struct recursive_segment_tree : protected recursive_range_query<segment_tree,O, Rank>
    {
        using child_type = recursive_range_query_operator<segment_tree, O, Rank - 1>;
        using base_type = segment_tree<child_type>;
        using value_type = O::value_type;
        using key_type = std::array<int,Rank>;
        using key_view_type = std::span<const int,Rank>;
        using tensor_type = std::vector<typename child_type::tensor_type>;

        friend class recursive_range_query_operator<segment_tree,O, Rank>;
    public:
        using base_type::base_type;

        recursive_segment_tree() = default;

        recursive_segment_tree(const tensor_type &elements) :
            base_type(std::vector(std::from_range,
                make_vector_optional<typename child_type::tensor_type,child_type>(elements)))
        {
        }

        template<typename RangeU, typename RangeV>
        value_type query(RangeU U, RangeV V)
        {
            key_view_type A(U),B(V);
            return query(A,B);
        }

        value_type query(key_view_type U,key_view_type V)
        {
            auto U_=U.template subspan<1,Rank-1>();
            auto V_=V.template subspan<1,Rank-1>();
            auto u = U[0], v=V[0];
            std::optional result(base_type::query(u,v));
            if (result.has_value()) return result.value().query(U_,V_);
            return O::neutral;
        }
    };

}

#endif //CPLIBRARY_RECURSIVE_H
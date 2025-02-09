//
// Created by ramizouari on 21/12/24.
//

#ifndef CPLIBRARY_PRODUCT_GROUP_H
#define CPLIBRARY_PRODUCT_GROUP_H
#include <tuple>

#include "algebra/structures.h"
namespace cp::groups
{
    namespace detail {
        template <typename Op, typename  Tp, std::size_t ... Is>
        auto opH2 (Op const & op, Tp const & t1, Tp const & t2,std::index_sequence<Is...> const &)
        {
            return std::make_tuple( op(std::get<Is>(t1), std::get<Is>(t2))... );
        }

        template <typename Op, typename Tp>
        auto opH1 (Op const & op, Tp const & t1, Tp const & t2)
        { return opH2(op, t1, t2,
                      std::make_index_sequence<std::tuple_size<Tp>{}>{});
        }
    }



    template<typename ...Ts>
    struct algebraic_product {
        std::tuple<Ts...> elements{};
        algebraic_product() = default;
        algebraic_product(Ts const&... elements) : elements(elements...) {}
        algebraic_product(std::tuple<Ts...> &&elements) : elements(std::move(elements)) {}
    };

    template <additive_monoid ... Ts>
    auto operator+ (algebraic_product<Ts...> const & t1, algebraic_product<Ts...> const & t2)
    { return detail::opH1(std::plus{}, t1.elements, t2.elements); }

    template <additive_group ... Ts>
    auto operator- (algebraic_product<Ts...> const & t1, algebraic_product<Ts...> const & t2)
    { return detail::opH1(std::minus{}, t1.elements, t2.elements); }

    template <multiplicative_monoid ... Ts>
    auto operator* (algebraic_product<Ts...> const & t1, algebraic_product<Ts...> const & t2)
    { return detail::opH1(std::multiplies{}, t1.elements, t2.elements); }

    template <multiplicative_group ... Ts>
    auto operator/ (algebraic_product<Ts...> const & t1, algebraic_product<Ts...> const & t2)
    { return detail::opH1(std::divides{}, t1.elements, t2.elements); }

    template <additive_monoid ... Ts>
    auto operator+= (algebraic_product<Ts...>  & t1, algebraic_product<Ts...> const & t2)
    { return t1=detail::opH1(std::plus{}, t1.elements, t2.elements); }

    template <additive_group ... Ts>
    auto operator-= (algebraic_product<Ts...> & t1, algebraic_product<Ts...> const & t2)
    { return t1=detail::opH1(std::minus{}, t1.elements, t2.elements); }

    template <multiplicative_monoid ... Ts>
    auto operator*= (algebraic_product<Ts...> & t1, algebraic_product<Ts...> const & t2)
    { return t1=detail::opH1(std::multiplies{}, t1.elements, t2.elements); }

    template <multiplicative_group ... Ts>
    auto operator/= (algebraic_product<Ts...> & t1, algebraic_product<Ts...> const & t2)
    { return t1=detail::opH1(std::divides{}, t1.elements, t2.elements); }

}

#endif //CPLIBRARY_PRODUCT_GROUP_H

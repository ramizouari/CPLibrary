//
// Created by ramizouari on 13/12/24.
//

#ifndef ALGEBRAIC_STRUCTURES_H
#define ALGEBRAIC_STRUCTURES_H
#include <concepts>
#include "types.h"

namespace cp
{
    template<typename G>
    concept additive_monoid = requires (G a,G b)
    {
      { a + b } -> std::convertible_to<G>;
      { a+=b } -> std::convertible_to<G>;
    };

    template<typename G>
    concept monoid = additive_monoid<G>;


    template<typename G>
    concept additive_group =  additive_monoid<G> && requires (G a,G b)
    {
        { a - b } -> std::convertible_to<G>;
        { a-=b } -> std::convertible_to<G>;
        { -a } -> std::convertible_to<G>;
    };

    template<typename G>
    concept group = additive_group<G>;

    template<typename G>
    concept abelian_group =  additive_group<G> && std::convertible_to<integer,G> && requires (G a,integer k)
    {
      { a + k } -> std::convertible_to<G>;
      { k + a } -> std::convertible_to<G>;
      { a +=k } -> std::convertible_to<G>;
      { a=k } -> std::convertible_to<G>;
      { a-k } -> std::convertible_to<G>;
      { k-a } -> std::convertible_to<G>;
      { a-=k } -> std::convertible_to<G>;
    };

    template<typename M,typename R>
    concept module = abelian_group<M> && requires (M a,R k)
    {
        { k*a } -> std::convertible_to<M>;
        { a*k } -> std::convertible_to<M>;
        { a*= k } -> std::convertible_to<M>;
    };

    template<typename M,typename R>
    concept vector_space = module<M,R> && requires (M a,R k)
    {
        { a/k } -> std::convertible_to<M>;
        { a/= k } -> std::convertible_to<M>;
    };


    template<typename G>
    concept multiplicative_monoid = requires (G a,G b)
    {
        { a * b } -> std::convertible_to<G>;
        { a*=b } -> std::convertible_to<G>;
    };

    template<typename G>
    concept multiplicative_group =  multiplicative_monoid<G> && requires (G a,G b)
    {
        { a / b } -> std::convertible_to<G>;
        { a/=b } -> std::convertible_to<G>;
    };


    template<typename G>
    concept abelian_multiplicative_group =  multiplicative_group<G> && std::convertible_to<integer,G> && requires (G a,integer k)
    {
        { a * k } -> std::convertible_to<G>;
        { k * a } -> std::convertible_to<G>;
        { a *=k } -> std::convertible_to<G>;
        { a=k } -> std::convertible_to<G>;
        { a/k } -> std::convertible_to<G>;
        { k/a } -> std::convertible_to<G>;
        { a/=k } -> std::convertible_to<G>;
    };

    template<typename R>
    concept ring = multiplicative_monoid<R> && group<R>;

    template<typename S>
    concept poset = requires (S a,S b)
    {
        { a <=> b } -> std::convertible_to<std::partial_ordering>;
    };

    template<typename S>
    concept ordered_set =  poset<S> && requires (S a,S b)
    {
        { a<=> b } -> std::convertible_to<std::weak_ordering>;
    };

    template<typename S>
    concept hashable_set =  requires (S a,std::hash<S> H)
    {
        {H(a)} -> std::convertible_to<std::uint64_t>;
    };


}

#endif //ALGEBRAIC_STRUCTURES_H

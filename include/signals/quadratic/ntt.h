//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_NTT_H
#define CPLIBRARY_NTT_H

#include <memory>
#include "rings/quadratic/dynamic.h"
#include "nt/modular_arithmetic.h"
#include "signals/fft.h"
#include "nt/number_theory.h"

namespace cp::signals
{
    template<integer m,natural o>
    struct cyclic_extension_t
    {
        using type=quadratic::dynamic_quadratic_extension<typename cyclic_extension_t<m,o-1>::type>;
        struct iterator
        {
            type holder;
        };
    };

    template<integer m>
    struct cyclic_extension_t<m,0>
    {
        using type=cyclic<m>;
    };
    template<natural m, natural o>
    using cyclic_extension=cyclic_extension_t<m,o>::type;


    template<integer n,quadratic::extension_type type,typename U,typename V>
    struct abstract_quadratic_ntt : public abstract_fft<quadratic::quadratic_extension_t<cyclic<n>,type,U,V>>, protected default_factoriser_t
    {
        using R=quadratic::quadratic_extension_t<cyclic<n>,type,U,V>;
        using abstract_fft<R>::transform;
        abstract_quadratic_ntt(std::shared_ptr<abstract_factoriser> _F=default_factoriser):F(_F){}
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
        virtual ~abstract_quadratic_ntt()= default;
    };
}

#endif //CPLIBRARY_NTT_H

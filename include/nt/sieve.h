//
// Created by ramizouari on 15/09/24.
//

#ifndef CPLIBRARY_SIEVE_H
#define CPLIBRARY_SIEVE_H
#include "number_theory.h"
#include "algebra/binary_operation.h"
#include <concepts>
#include <utility>

namespace cp
{

    template<typename F, typename T>
    concept arithmetic_generator=  requires(F f, int i)
    {
        {f(i) } -> std::convertible_to<T>; // Requires f(i) to return a T
    };

    template<typename F, typename T>
    concept full_arithmetic_generator=  requires(F f, int i, const abstract_factoriser & factoriser)
    {
        {f(i,factoriser) } -> std::convertible_to<T>; // Requires f(i) to return a T
    };

    template<typename F, typename T>
    concept grouped_arithmetic_generator=  requires(F f, std::vector<int> I)
    {
        {f(I) } -> std::convertible_to<std::vector<T>>; // Requires f(i) to return a T
    };

    template<typename F, typename T>
    concept grouped_full_arithmetic_generator=  requires(F f, std::vector<int> I, const abstract_factoriser & factoriser)
    {
        {f(I,factoriser) } -> std::convertible_to<std::vector<T>>; // Requires f(i) to return a T
    };


    enum sieve_mode : bool
    {
        Primes,
        PrimePowers
    };

    enum sieve_order
    {
        Natural,
        Prime,
        Multiplicity
    };

    template<typename I>
    void sieve_inplace(std::vector<I> & S, const cp::abstract_factoriser& F,
                       const binary_operation<I> & Op = multiplies_t<I>{},
                       sieve_mode mode = sieve_mode::PrimePowers)
    {
        for(int k=2;k<S.size();k++)
        {
            cp::integer q;
            if(mode == sieve_mode::Primes)
                q= F.smallest_divisor(k);
            else
                q = F.smallest_maximal_prime_power_divisor(k);
            S[k] = Op(S[q],S[k/q]);
        }
    }

    template<typename I,arithmetic_generator<I> Generator>
    void sieve_inplace(std::vector<I> & S, Generator && generator, const cp::abstract_factoriser& F,
                       const binary_operation<I> & Op = multiplies_t<I>{},
                       sieve_mode mode = sieve_mode::PrimePowers)
    {
        for(int k=2;k<S.size();k++)
        {
            cp::integer q;
            if(mode == sieve_mode::Primes)
                q= F.smallest_divisor(k);
            else
                q = F.smallest_maximal_prime_power_divisor(k);
            if(q==k)
                S[k] = generator(k);
            else S[k] = Op(S[q],S[k/q]);
        }
    }

    template<typename I,full_arithmetic_generator<I> Generator>
    void sieve_inplace(std::vector<I> & S, Generator && generator, const cp::abstract_factoriser& F,
                       const binary_operation<I> & Op = multiplies_t<I>{},
                       sieve_mode mode = sieve_mode::PrimePowers)
    {
        for(int k=2;k<S.size();k++)
        {
            cp::integer q;
            if(mode == sieve_mode::Primes)
                q= F.smallest_divisor(k);
            else
                q = F.smallest_maximal_prime_power_divisor(k);
            if(q==k)
                S[k] = generator(k,F);
            else S[k] = Op(S[q],S[k/q]);
        }
    }

    template<typename I>
    struct registered_sieve_array
    {
        std::span<I> content;
        binary_operation_ptr<I> op;
        registered_sieve_array(std::span<I> content, binary_operation_ptr<I> op): content(content), op(op){}
    };


    template<typename I = integer>
    class sieve
    {
        std::array<std::vector<registered_sieve_array<I>>,2> registered;
        const abstract_factoriser & F;
        std::array<std::vector<integer>,2> factors_list;
        integer n;
    public:
        explicit sieve(const abstract_factoriser & F,integer n,sieve_order order = Natural) :F(F), n(n)
        {
            for(int i=2;i<=n;i++)
            {
                if (i==F.smallest_divisor(i)) factors_list[Primes].push_back(i);
                if(i==F.smallest_maximal_prime_power_divisor(i)) factors_list[PrimePowers].push_back(i);
            }
            switch(order)
            {
                case Prime:
                    for(auto& factors: factors_list)
                        std::sort(factors.begin(),factors.end(),[&F](auto x,auto y){
                            return std::make_pair(F.smallest_divisor(x),x) < std::make_pair(F.smallest_divisor(y),y);
                        });
                    break;
                case Multiplicity:
                    for(auto& factors: factors_list)
                        std::sort(factors.begin(),factors.end(),[&F](auto x,auto y){
                            return std::make_pair(F.smallest_prime_multiplicity(x),x) < std::make_pair(F.smallest_prime_multiplicity(y),y);
                        });
                    break;
                default:
                    break;
            }
        }

        template<grouped_arithmetic_generator<I> Generator>
        void add(std::span<I> S, Generator && generator,
                 std::shared_ptr<binary_operation<I>> op  = std::make_shared<multiplies_t<I>>())
        {
            S[1] = op->neutral_element();
            for(auto p:factors_list[Primes])
            {
                std::vector<int> J;
                for(integer q = 1; q <= n;q*=p)
                    J.push_back(q);
                auto H=generator(J);
                for(int i=1;i<J.size();i++)
                    S[J[i]] = H[i];
            }
            registered[PrimePowers].emplace_back(S,op);
        }

        template<grouped_full_arithmetic_generator<I> Generator>
        void add(std::span<I> S, Generator && generator,
                 std::shared_ptr<binary_operation<I>> op  = std::make_shared<multiplies_t<I>>())
        {
            S[1] = op->neutral_element();
            for(auto p:factors_list[Primes])
            {
                std::vector<int> J;
                for(integer q = 1; q <= n;q*=p)
                    J.push_back(q);
                auto H=generator(J,F);
                for(int i=1;i<J.size();i++)
                    S[J[i]] = H[i];
            }
            registered[PrimePowers].emplace_back(S,op);
        }


        template<arithmetic_generator<I> Generator>
        void add(std::span<I> S, Generator && generator,
                 std::shared_ptr<binary_operation<I>> op  = std::make_shared<multiplies_t<I>>(),
                 sieve_mode mode = sieve_mode::PrimePowers)
        {
            S[1] = op->neutral_element();
            for(auto q:factors_list[mode])
                S[q] = generator(q);
            registered[mode].emplace_back(S,op);
        }

        template<full_arithmetic_generator<I> Generator> requires (!grouped_full_arithmetic_generator<Generator,I>)
        void add(std::span<I> S, Generator && generator,
                 std::shared_ptr<binary_operation<I>> op  = std::make_shared<multiplies_t<I>>(),
                 sieve_mode mode = sieve_mode::PrimePowers)
        {
            S[1] = op->neutral_element();
            for(auto q:factors_list[mode])
                S[q] = generator(q,F);
            registered[mode].emplace_back(S,op);
        }

        void generate()
        {
            for(int i=2;i<=n;i++)
            {
                auto d = F.smallest_divisor(i);
                auto q = F.smallest_maximal_prime_power_divisor(i);
                for(auto & [S,op]: registered[Primes]) if(i!=q)
                    S[i] = op(S[d],S[i/d]);
                for(auto & [S,op]: registered[PrimePowers]) if(i!=q)
                    S[i] = op(S[q],S[i/q]);
            }
        }

    };

}

#endif //CPLIBRARY_SIEVE_H

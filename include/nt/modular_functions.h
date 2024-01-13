//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_MODULAR_FUNCTIONS_H
#define CPLIBRARY_MODULAR_FUNCTIONS_H
#include "modular_arithmetic.h"
#include "polynomial/ring_extension.h"
#include "functions.h"
#include <random>

template<cp::integer m>
struct std::hash<cp::cyclic<m>>
{
    inline static std::random_device dev=std::random_device();
    inline static std::mt19937 g=std::mt19937(dev());
    inline static constexpr cp::integer M=1e9+7;
    std::uniform_int_distribution<cp::integer> d=std::uniform_int_distribution<cp::integer>(1,M);
    cp::integer a=d(g),b=d(g);
    public:
    size_t operator()(const cp::cyclic<m> &x) const noexcept
    {
        return (a*static_cast<cp::integer>(x)+b)%M;
    }
};

namespace cp
{
    template<typename cyclic_ring>
    integer discrete_log(cyclic_ring a, cyclic_ring r)
    {
        integer s=std::ceil(std::sqrt(cyclic_ring::m));
        cyclic_ring u=pow(a,s),w=1;
        std::unordered_map<cyclic_ring,integer> mapper;
        for(integer i=0;i<=s;i++,w*=a)
            mapper[r*w]=i;
        w=u;
        for(integer i=1;i<=s;i++,w*=u)
            if(mapper.count(w))
                return i*s-mapper[w];
        return -1;
    }

    inline std::vector<integer> inverse_table(int n,int prime)
    {
        std::vector<integer> I(n + 1);
        I[0] = I[1] = 1;
        for (int i = 2; i <= n; i++)
            I[i] = I[prime % i] *
                   (prime - prime / i) % prime;
        return I;
    }

    inline integer primitive_root_of_unity(integer p,abstract_factoriser &F)
    {
        auto phi=carmichael_totient(p,F);
        auto D=F.divisors_list(phi);
        for(integer k=2;k<p-1;k++) if(std::gcd(k,p)==1)
        {
            bool is_primitive=true;
            for (auto d: D)
                if(d< phi && pow(cyclic<dynamic_modulus>(k,p),d,p)==1)
                {
                    is_primitive=false;
                    break;
                }
            if(is_primitive)
                return k;
        }
        return 0;
    }

    template <integer p>
    integer primitive_root_of_unity(abstract_factoriser& F)
    {
        static auto phi = totient(p,F);
        static auto D = F.divisors_list(phi);
        for (integer k = 2; k < p - 1; k++)
        {
            bool is_primitive = true;
            for (auto d : D)
                if (d < phi && pow<d_cyclic>(k, d) == 1)
                {
                    is_primitive = false;
                    break;
                }
            if (is_primitive)
                return k;
        }
        return 0;
    }

    template<integer m>
    integer legendre_symbol(cyclic<m> a)
    {
        integer r;
        if constexpr (m==dynamic_modulus)
            r= (integer) pow(a, (a.modulus() - 1) / 2,a.modulus());
        else
            r= (integer)pow(a, (a.modulus() - 1) / 2);
        if (r > a.modulus() / 2)
            r -= a.modulus();
        return r;
    }

    template<integer m>
    cyclic<m> sqrt(cyclic<m> n)
    {
        using cyclic_field = cyclic<m>;
        if (n.modulus() == 2)
            return n;
        cyclic_field a = 2;
        while (legendre_symbol(a*a-n) != -1)
            ++a;
        extension_polynomial_t<cyclic_field> q = { polynomial<cyclic_field>({n-a * a,0,1}) };
        d_ring_extension<cyclic_field> phi(std::vector<cyclic_field>{ a,1 }, q);
        return pow(phi, (cyclic_field::m+1)/2,q)[0];
    }
}

#endif //CPLIBRARY_MODULAR_FUNCTIONS_H

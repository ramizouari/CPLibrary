//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_QUADRATIC_EXTENSION_FIXED_H
#define CPLIBRARY_QUADRATIC_EXTENSION_FIXED_H
#include "algebra/abstract_algebra.h"

namespace cp::quadratic
{
    enum class extension_type
    {
        dynamic_extension,
        static_extension,
        fixed_extension
    };
    template<typename R,extension_type type,typename A1,typename A2>
    struct quadratic_extension_t;

    template<typename R,typename A1,A1 a, typename A2, A2 b>
    struct quadratic_extension_t<R, extension_type::fixed_extension,std::integral_constant<A1,a>,std::integral_constant<A2,b>>
    {
        std::array<R,2> A;
        R& operator[](int i)
        {
            return A[i];
        }
        template<typename R1,typename R2>
        quadratic_extension_t(R1 x,R2 y):A{x,y}{}
        quadratic_extension_t(integer x=0):A{x,R{}}{}
        quadratic_extension_t(R x, R y):A{x,y}{}
        quadratic_extension_t operator+(const quadratic_extension_t &x) const
        {
            return {A[0]+x.A[0],A[1]+x.A[1]};
        }
        quadratic_extension_t operator-(const quadratic_extension_t &x) const
        {
            return {A[0]-x.A[0],A[1]-x.A[1]};
        }
        quadratic_extension_t operator*(const quadratic_extension_t &x) const
        {
            return {b*A[1]*x.A[1]+A[0]*x.A[0],A[0]*x.A[1]+A[1]*x.A[0]+a*A[1]*x.A[1]};
        }

        quadratic_extension_t conj() const
        {
            return {A[0]-a*A[0],-A[1]};
        }

        quadratic_extension_t norm() const
        {
            return A[0]*A[0]+a*A[0]*A[1]-b*A[1]*A[1];
        }

        quadratic_extension_t inv() const
        {
            auto n=norm();
            return conj()/n;
        }

        bool operator==(const quadratic_extension_t &x) const = default;

        quadratic_extension_t& operator+=(const quadratic_extension_t &x)
        {
            A[0]+=x.A[0];
            A[1]+=x.A[1];
            return *this;
        }

        quadratic_extension_t& operator-=(const quadratic_extension_t &x)
        {
            A[0]-=x.A[0];
            A[1]-=x.A[1];
            return *this;
        }

        quadratic_extension_t& operator*=(const quadratic_extension_t &x)
        {
            return *this=(*this)*x;
        }

        quadratic_extension_t& operator/=(const quadratic_extension_t &x)
        {
            return *this=(*this)*x.inv();
        }

        quadratic_extension_t operator-() const
        {
            return {-A[0],-A[1]};
        }

        quadratic_extension_t operator/(const quadratic_extension_t &x) const
        {
            return (*this)*x.inv();
        }

        quadratic_extension_t operator~() const
        {
            return conj();
        }
        static std::array<R,2> modulus()
        {
            return {R(a),R(b)};
        }
    };
    template<typename R, integer a,integer b>
    using integral_quadratic_extension=quadratic_extension_t<R,extension_type::fixed_extension,std::integral_constant<integer,a>,std::integral_constant<integer,b>>;
    template<typename R,R a,R b>
    using quadratic_extension=quadratic_extension_t<R,extension_type::fixed_extension,std::integral_constant<R,a>,std::integral_constant<R,b>>;
}

#endif //CPLIBRARY_FIXED_H

//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_QUADRATIC_STATIC_H
#define CPLIBRARY_QUADRATIC_STATIC_H
#include "fixed.h"
namespace cp::quadratic
{
    template<typename R,extension_type type,typename A1,typename A2>
    struct quadratic_extension_t;

    template<typename R>
    struct quadratic_extension_t<R, extension_type::static_extension, void, void>
    {
        std::array<R,2> A;
        R& operator[](int i)
        {
            return A[i];
        }
        template<typename R1,typename R2>
        quadratic_extension_t(R1 x,R2 y):A{x,y}{}
        template<std::integral I>
        quadratic_extension_t(I x): A{x,R{}}{}
        quadratic_extension_t() : A{} {}
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

        R norm() const
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

        quadratic_extension_t operator/(const R &x) const
        {
            return {A[0]/x,A[1]/x};
        }


        quadratic_extension_t operator~() const
        {
            return conj();
        }
        static std::array<R,2> modulus()
        {
            return {a,b};
        }

        static void set_modulus(const std::array<R,2> &M_)
        {
            a=M_[0];
            b=M_[1];
        }

        static void set_modulus(R a_,R b_)
        {
            a=a_;
            b=b_;
        }

        inline static R a,b;
    };

    template<typename R>
    using static_quadratic_extension=quadratic_extension_t<R,extension_type::static_extension,void,void>;
}

#endif //CPLIBRARY_STATIC_H

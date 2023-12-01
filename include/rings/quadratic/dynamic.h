//
// Created by ramizouari on 29/11/23.
//

#ifndef CPLIBRARY_QUADRATIC_EXTENSION_DYNAMIC_H
#define CPLIBRARY_QUADRATIC_EXTENSION_DYNAMIC_H
#include "fixed.h"
namespace cp::quadratic
{
    template<typename R,extension_type type,typename A1,typename A2>
    struct quadratic_extension_t;

    template<typename R>
    struct quadratic_extension_t<R, extension_type::dynamic_extension, void, void>
    {
        std::array<R,2> A,B{};
        bool definite=false;
        R& operator[](int i)
        {
            return A[i];
        }
        template<typename R1,typename R2>
        quadratic_extension_t(R1 x,R2 y,std::array<R,2> B_):A{x,y},B(B_),definite(true){}
        quadratic_extension_t(integer x=0):A{x,R{}} {}
        quadratic_extension_t(R x, R y,std::array<R,2> B_):A{x,y},B(B_),definite(true){}
        quadratic_extension_t operator+(const quadratic_extension_t &x) const
        {
            return {A[0]+x.A[0],A[1]+x.A[1],B};
        }
        quadratic_extension_t operator-(const quadratic_extension_t &x) const
        {
            return {A[0]-x.A[0],A[1]-x.A[1],B};
        }
        quadratic_extension_t operator*(const quadratic_extension_t &x) const
        {
            auto &B=definite?this->B:x.B;
            return {B[1]*A[1]*x.A[1]+A[0]*x.A[0],A[0]*x.A[1]+A[1]*x.A[0]+B[0]*A[1]*x.A[1],B};
        }

        quadratic_extension_t conj() const
        {
            return {A[0]-B[0]*A[0],-A[1],B};
        }

        quadratic_extension_t norm() const
        {
            return A[0]*A[0]+B[0]*A[0]*A[1]-B[1]*A[1]*A[1];
        }

        quadratic_extension_t inv() const
        {
            auto n=norm();
            return conj()/n;
        }

        bool operator==(const quadratic_extension_t &x) const = default;

        quadratic_extension_t& operator+=(const quadratic_extension_t &x)
        {
            if(!definite) {
                B = x.B;
                definite=x.definite;
            }
            A[0]+=x.A[0];
            A[1]+=x.A[1];
            return *this;
        }

        quadratic_extension_t& operator-=(const quadratic_extension_t &x)
        {
            if(!definite) {
                B = x.B;
                definite=x.definite;
            }
            A[0]-=x.A[0];
            A[1]-=x.A[1];
            return *this;
        }

        quadratic_extension_t& operator*=(const quadratic_extension_t &x)
        {
            if(!definite) {
                B = x.B;
                definite=x.definite;
            }
            return *this=(*this)*x;
        }

        quadratic_extension_t& operator/=(const quadratic_extension_t &x)
        {
            if(!definite) {
                B = x.B;
                definite=x.definite;
            }
            return *this=(*this)*x.inv();
        }

        quadratic_extension_t operator-() const
        {
            return {-A[0],-A[1],B};
        }

        quadratic_extension_t operator/(const quadratic_extension_t &x) const
        {
            return (*this)*x.inv();
        }

        quadratic_extension_t operator~() const
        {
            return conj();
        }
        std::array<R,2> modulus()
        {
            return B;
        }

        void set_modulus(const std::array<R,2> &M_)
        {
            B[0]=M_[0];
            B[1]=M_[1];
        }

        void set_modulus(R a_,R b_)
        {
            B[0]=a_;
            B[1]=b_;
        }
    };

    template<typename R>
    using dynamic_quadratic_extension=quadratic_extension_t<R,extension_type::dynamic_extension,void,void>;
}

#endif //CPLIBRARY_FIXED_H

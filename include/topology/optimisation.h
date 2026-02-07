//
// Created by ASUS on 01/12/2021.
//

#ifndef __OPTIMISATION__
#define __OPTIMISATION__
#include "linalg/matrix.h"
#include "../algebra/abstract_algebra.h"
#include "topology.h"
#include <functional>
namespace cp::topology
{
    using namespace linalg;

    template<typename E,real_type IK=IR,typename Norm=L2_inner_product<real,E>>
    class gradient_descent
    {
        inline static constexpr Norm N{};
    protected:
        IK p=.05;
        IK eps;
        derivator<E,real,E>& D;
    public:
        gradient_descent(derivator<E,IK,E> &d,IK _eps=1e-3):D(d),eps(_eps) {}
        E argmin(const std::function<IK(E)>& f,E x) const
        {
            for (; N.norm(D.gradient(f, x)) > eps; x -= p * D.gradient(f, x));
            return x;
        }
        E argmin(const std::function<IK(E)>& f, E x,int L) const
        {
            for (; N.norm(D.gradient(f, x)) > eps && L--; x -= p * D.gradient(f, x));
            return x;
        }
        E argmin(const std::function<std::pair<IK,E>(E)>& f, E x) const
        {
            auto P = f(x);
            while (N.norm(P.second) > eps)
            {
                x -= p * P.second;
                P = f(x);
            }
            return x;
        }

        E argmin(const std::function<std::pair<real, E>(E)>& f, E x,int L) const
        {
            auto P = f(x);
            while (N.norm(P.second) > eps && L--)
            {
                x -= p * P.second;
                P = f(x);
            }
            return x;
        }

    };

    template<typename E,real_type IK=IR,typename InnerProduct=L2_inner_product<IK,E>>
    class barzilai_borwein_gradient_descent
    {
        IK p=.1;
        IK eps=1e-8;
        derivator<E,real,E>& D;
        inline static constexpr InnerProduct B{};
    public:
        barzilai_borwein_gradient_descent(derivator<E, IK,E>& d, real _p):D(d),p(_p){}

        E argmin(const std::function<IK(E)>& f, E s,int L)
        {
            this->p = 0.1;
            E x = s- this->p*this->D.gradient(f, s);
            for (; B.norm(this->D.gradient(f, x)) > this->eps && L; x -= this->p * this->D.gradient(f, x))
            {
                update_rate(f, x,s);
                s = x;
                L--;
            }
            return x;
        }

        E argmin(const std::function<std::pair<IK, E>(E)>& f, E s, int L)
        {
            auto P = f(s);
            E x = s - this->p * f(s).second;
            auto Q = f(x);
            while (B.norm(P.second) > eps && L--)
            {
                update_rate(x, s,Q.second,P.second);
                s = x;
                x -= this->p * Q.second;
                P = std::move(Q);
                Q = f(x);
            }
            return x;
        }


        virtual void update_rate(const std::function<IK(E)>& f, const E& x,const E& s)
        {
            auto L = this->D.gradient(f, x) - this->D.gradient(f, s);
            this->p = B.inner_product(L,x - s) / B.inner_product(L,L);
        }

        virtual void update_rate(const E& x, const E& s,const E& dx,const E& ds)
        {
            auto L = dx - ds;
            this->p = B.inner_product(L, x - s) / B.inner_product(L, L);
        }
    };
}
#endif

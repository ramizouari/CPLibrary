//
// Created by ramizouari on 01/12/2021.
//
#ifndef __POLYNOMIAL__H__
#define __POLYNOMIAL__H__
#include <vector>
#include <map>
#include "algebra/abstract_algebra.h"

namespace cp
{
    /**
* @brief Polynomial class
* @details This is the class of polynomials over a commutative ring R
* @tparam R the type of the coefficients
* @Requirements
* <strong>R</strong> is a commutative ring
* @Notes
* Formally this class is simply R[x]
*/
    template<typename R>
    struct polynomial
    {
        std::vector<R> p;

    public:

        const std::vector<R>& data() const
        {
            return p;
        }
        std::vector<R>& data()
        {
            return p;
        }

        bool operator==(const R& a) const
        {
            if (is_zero(a))
                return p.empty();
            else return degree() == 0 && p.front() == a;
        }

        bool operator==(const polynomial<R> &) const = default;

        polynomial& reduce()
        {
            while(!p.empty() && is_zero(p.back()))
                p.pop_back();
            return *this;
        }
        polynomial(R k):p(1,k)
        {
            reduce();
        }

        template<std::integral T>
        polynomial(T k):p(1,k)
        {
            reduce();
        }

        polynomial()=default;

        polynomial(std::vector<R> _p):p(std::move(_p))
        {
            reduce();
        }

        int degree() const
        {
            return p.size()-1;
        }

        int dim() const
        {
            return p.size();
        }

        auto& operator+=(const polynomial &O)
        {
            int r=std::max(p.size(),O.p.size());
            p.resize(r);
            for(int i=0;i<O.p.size();i++)
                p[i]+=O.p[i];
            reduce();
            return *this;
        }

        auto& operator-=(const polynomial &O)
        {
            int r=std::max(p.size(),O.p.size());
            p.resize(r);
            for(int i=0;i<O.p.size();i++)
                p[i]-=O.p[i];
            reduce();
            return *this;
        }

        polynomial operator*(const polynomial &O) const
        {
            if(O.p.empty() || p.empty())
                return polynomial(0);
            int n=degree(),m=O.degree();
            polynomial q;
            q.p.resize(n+m+1);
            for(int i=0;i<=n;i++) for(int j=0;j<=m;j++)
                    q.p[i+j]+=p[i]*O.p[j];
            q.reduce();
            return q;
        }

        auto& operator*=(const polynomial &O)
        {
            auto r=(*this)*O;
            p.swap(r.p);
            return *this;
        }

        auto operator+(const polynomial &O) const
        {
            auto r=*this;
            return r+=O;
        }

        auto operator-(const polynomial &O) const
        {
            auto r=*this;
            return r-=O;
        }

        auto operator-() const
        {
            auto r=*this;
            for(auto &s:r.p)
                s=-s;
            return r;
        }

        auto operator*=(R a)
        {
            if(is_zero(a))
                p.clear();
            else for(auto& s:p)
                    s*=a;
            reduce();
            return *this;
        }

        bool operator!=(R a) const
        {
            return !(*this==a);
        }

        auto& operator+=(R a)
        {
            return *this+=polynomial({a});
        }

        auto& operator-=(R a)
        {
            return *this-=polynomial({a});
        }

        auto operator+(R a) const
        {
            auto q=*this;
            return q+=a;
        }

        auto operator-(R a) const
        {
            auto q=*this;
            return q-=a;
        }

        /**
        * @details creates a preorder between polynomials based on the degree
        * @Requirements:
        * None
        * @Notes
        * This function is essential for the euclidean algorithm to work
        */
        bool operator<(const polynomial &O) const
        {
            return degree() < O.degree();
        }

        /**
         * @brief Polynomial self division
        * @details Divides the polynomial by a constant and stores the result in itself
        * @Requirements
        * One of the following:
         * <ul>
        * <li> R is an integral ring [2]
        * <li> k is invertible
         * </ul>
        */

        auto& operator/=(R k)
        {
            for(auto &s:p)
                s/=k;
            return *this;
        }

        auto operator/(R k) const
        {
            auto q=*this;
            return q/=k;
        }

        /**
         * @brief Euclidean division
        * @details Applies euclidean division between two polynomials
        * @Requirements
        * One of the following
        * <ul>
        * <li> R is a field [1]
        * <li> R is an euclidean domain [2]
        * <li> R is a commutative ring, and the dominant coefficient of <strong>O</strong> is inversible
        * <ul/>
        * @Notes
        * Even that condition [1] is a special case of [2], given that some properties of euclidean division are
        * guaranteed only if <strong>R</strong> is a field, We will seperate the two cases
        */
        std::pair<polynomial,polynomial> euclidean_division(const polynomial &O) const
        {
            if(degree() < O.degree())
                return {R(0),*this};
            polynomial q,r=*this;
            int n=degree(),m=O.degree();
            q.p.resize(n-m+1);
            for(int i=n;i>=m;i--)
            {
                q.p[i-m]=r[i]/O.p[m];
                for(int j=0;j<=m;j++)
                    r.p[i+j-m]-=q.p[i-m]*O.p[j];
            }
            r.reduce();
            return {q,r};
        }

        //Get the quotient of the euclidean division
        polynomial operator/(const polynomial &O) const
        {
            return euclidean_division(O).first;
        }

        //Get the remainder of the euclidean division
        polynomial operator%(const polynomial &O) const
        {
            return euclidean_division(O).second;
        }

        auto& operator/=(polynomial &O)
        {
            p.swap(((*this)/O).p);
            return *this;
        }

        auto &operator[](int k)
        {
            return p[k];
        }

        const auto& operator[](int k) const
        {
            return p[k];
        }

        polynomial derivative() const
        {
            if (p.empty())
                return {};
            polynomial D;
            D.p.resize(degree());
            for(int i=1;i<=degree();i++)
                D.p[i-1]=p[i]*R(i);
            return D;
        }

        /**
        * @brief Polynomial Evaluation
        * @details Evaluates the polynomial at a given point
        * Evaluates the polynomial over an associative R-algebra H
        * @Requirements
        * <strong>H</strong> is an associative algebra over <strong>R</strong>
        */
        template<typename H>
        std::common_type<H,R>::type operator()(H a) const
        {
            typename std::common_type<H,R>::type r{};
            for(int i=degree();i>=0;i--)
                r=r*a+p[i];
            return r;
        }

        auto begin()
        {
            return p.begin();
        }

        auto end()
        {
            return p.end();
        }

        auto begin() const
        {
            return p.begin();
        }

        auto end() const
        {
            return p.end();
        }

        explicit operator std::vector<R>&()
        {
            return p;
        }

        explicit operator const std::vector<R>&() const
        {
            return p;
        }
    };

    template<typename R>
    polynomial<R> operator*(R a,const polynomial<R> &p)
    {
        auto q=p;
        return q*=a;
    }

    template<typename R>
    polynomial<R> operator-(R a,const polynomial<R> &p)
    {
        auto q=-p;
        return q+=a;
    }

/**
* @brief The functional identity polynomial
* @details This constant is the generator of all polynomials over R.
* @Notes
* Formally, it is simply the polynomial X:x->x
*/
    template<typename R>
    const polynomial<R> X=polynomial<R>(std::vector<R>{0,1});

/**
* @brief Newton Interpolation
* @details Applies Lagrange Interpolation over points (x,y) using Newton's method
* @Requirements
 * <ol>
* <li> x does not have a duplicate element
* <li> One of the following
 * <ul>
*   <li> R is a field
*   <li> (s-t) is invertible for all elements s,t in x
 *  </ul>
 * </ol>
*/
    template<typename R>
    polynomial<R> newton_interpolation(const std::vector<R> &x,const std::vector<R> &y)
    {
        int n=x.size()-1;
        std::vector<std::vector<R>> d(n+1,std::vector<R>(n+1));
        for(int i=0;i<=n;i++)
            d[i][i]=y[i];
        for(int r=1;r<=n;r++) for(int i=0;i+r<=n;i++)
                d[i][i+r]=(d[i+1][i+r]-d[i][i+r-1])/(x[i+r]-x[i]);
        polynomial<R> p,u=R(1);
        for(int i=0;i<=n;i++) {
            p +=d[0][i]*u;
            u*=(X<R>-x[i]);
        }
        return p;
    }
}

#endif
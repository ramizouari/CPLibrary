//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_SPARSE_POLYNOMIAL_H
#define CPLIBRARY_SPARSE_POLYNOMIAL_H
#include <map>
#include <vector>
#include <cstdint>

namespace cp
{
    using integer=std::int64_t;
/**
 * @brief Sparse Polynomial
* @details This is the class of sparse polynomials over commutative ring R
* @Requirements
* <strong>R</strong> is a commutative ring
* @Recommendation
* <ol> <li> The coefficients are sparse. Formally a k-sparse polynomial p of degree n is a polynomial where:
* (card supp {p_1,..,p_n}) / n <= k
* <li> It is recommended that k<=0.01 </ol>
* @Notes
* Formally this class is simply R[x]
*/
    template<typename R>
    class sparse_polynomial
    {
        std::map<cp::integer,R> p;
        void reduce()
        {
            std::vector<cp::integer> to_del;
            for(auto [k,x]:p)
                if(is_zero(x))
                    to_del.push_back(k);
            for(auto k:to_del)
                p.erase(k);
        }
    public:
        sparse_polynomial(R k=R{})
        {
            p[0]=k;
            reduce();
        }

        template<std::integral I>
        sparse_polynomial(I k)
        {
            p[0]=k;
            reduce();
        }


        sparse_polynomial(const std::vector<R> &_p)
        {
            for(int i=0;i<_p.size();i++)
                p[i]=_p[i];
            reduce();
        }

        int degree() const
        {
            return p.empty()?-1:p.rbegin()->first;
        }

        sparse_polynomial& operator+=(const sparse_polynomial& O)
        {
            for(const auto& [k,s]:O.p)
            {
                p[k] += s;
                if(is_zero(p[k]))
                    p.erase(k);
            }
            return *this;
        }

        sparse_polynomial& operator-=(const sparse_polynomial& O)
        {
            for(const auto& [k,s]:O.p)
            {
                p[k] -= s;
                if(is_zero(p[k]))
                    p.erase(k);
            }
            return *this;
        }

        sparse_polynomial operator*(const sparse_polynomial &O) const
        {

            sparse_polynomial q;
            for(auto [i,u]:p) for(auto [j,v]:O.p)
                {
                    q.p[i+j]+=u*v;
                    if(is_zero(q.p[i+j]))
                        q.p.erase(i+j);
                }
            return q;
        }

        sparse_polynomial& operator*=(const sparse_polynomial &O)
        {
            auto r=(*this)*O;
            p.swap(r.p);
            return *this;
        }

        sparse_polynomial operator+(const sparse_polynomial &O) const
        {
            auto r=*this;
            return r+=O;
        }

        sparse_polynomial operator-(const sparse_polynomial &O) const
        {
            auto r=*this;
            return r-=O;
        }

        sparse_polynomial operator-() const
        {
            auto r=*this;
            for(auto &[_,s]:r.p)
                s=-s;
            return r;
        }

        sparse_polynomial operator*=(R a)
        {
            if(is_zero(a))
                p.clear();
            else for(auto& [_,s]:p)
                    s*=a;
            reduce();
            return *this;
        }

        sparse_polynomial& operator/=(R k)
        {
            for(auto &s:p)
                s/=k;
            return *this;
        }

        sparse_polynomial operator/(R k) const
        {
            auto q=*this;
            return q/=k;
        }

        R &operator[](cp::integer k)
        {
            return p[k];
        }

        const R& operator[](cp::integer k) const
        {
            return p.at(k);
        }

        /**
         * @brief Polynomial evaluation
        * @details Evaluates the polynomial at a point a.
        * @Requirements:
        * H is an associative algebra over R
        */
        template<typename H>
        std::common_type<H,R>::type operator()(H a) const
        {
            typename std::common_type<H,R>::type r=0,u=1;
            cp::integer i=0;
            for(auto [k,x]:p)
            {
                u*=pow(a,k-i);
                r+=u*x;
                i=k;
            }
            return r;
        }

        operator std::map<cp::integer, R>& ()
        {
            return p;
        }

        operator const std::map<cp::integer, R>& () const
        {
            return p;
        }

        auto size() const
        {
            return p.size();
        }

        std::map<cp::integer,R>& data()
        {
            return p;
        }

        const std::map<cp::integer,R>& data() const
        {
            return p;
        }

        sparse_polynomial derivative() const
        {
            sparse_polynomial q;
            for(auto [k,x]:p)
                q.p.emplace_hint(q.p.end(),k-1,k*x);
            return q;
        }
    };

    template<typename R>
    sparse_polynomial<R> Z=sparse_polynomial<R>(std::vector<R>{0,1});

    template<typename R>
    sparse_polynomial<R> operator*(R a,const sparse_polynomial<R> &p)
    {
        auto q=p;
        return q*=a;
    }

    template<typename R>
    sparse_polynomial<R> operator-(R a,const sparse_polynomial<R> &p)
    {
        auto q=-p;
        return q+=a;
    }

    template<typename R>
    sparse_polynomial<R> operator+(R a,const sparse_polynomial<R> &p)
    {
        auto q=p;
        return q+=a;
    }

}

#endif //CPLIBRARY_SPARSE_POLYNOMIAL_H
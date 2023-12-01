//
// Created by ramizouari on 28/11/23.
//

#ifndef CPLIBRARY_SPARSE_POLYNOMIAL_H
#define CPLIBRARY_SPARSE_POLYNOMIAL_H
#include <map>
#include <vector>
namespace cp
{

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
        std::map<int,R> p;
        void reduce()
        {
            std::vector<int> to_del;
            for(auto [k,x]:p)
                if(is_zero(x))
                    to_del.push_back(k);
            for(auto k:to_del)
                p.erase(k);
        }
    public:
        sparse_polynomial(R k=0)
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

        auto& operator+=(sparse_polynomial O)
        {
            for(const auto& [k,s]:O.p)
            {
                p[k] += O.p[k];
                if(is_zero(p[k]))
                    p.erase(k);
            }
            return *this;
        }

        auto& operator-=(sparse_polynomial O)
        {
            for(const auto& [k,s]:O.p)
            {
                p[k] -= O.p[k];
                if(is_zero(p[k]))
                    p.erase(k);
            }
            return *this;
        }

        auto operator*(const sparse_polynomial &O) const
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

        auto& operator*=(const sparse_polynomial &O)
        {
            auto r=(*this)*O;
            p.swap(r.p);
            return *this;
        }

        auto operator+(const sparse_polynomial &O) const
        {
            auto r=*this;
            return r+=O;
        }

        auto operator-(const sparse_polynomial &O) const
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

        auto& operator+=(R a)
        {
            return *this+=sparse_polynomial({a});
        }

        auto& operator-=(R a)
        {
            return *this+=sparse_polynomial({a});
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

        auto &operator[](int k)
        {
            return p[k];
        }

        const auto& operator[](int k) const
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
        H operator()(H a)
        {
            H r=0,u=1,i=0;
            for(auto [k,x]:p)
            {
                u*=pow(a,k-i);
                r+=u*x;
                i=k;
            }
        }

        operator std::map<int, R>& ()
        {
            return p;
        }

        operator const std::map<int, R>& () const
        {
            return p;
        }

        auto size() const
        {
            return p.size();
        }
    };

    template<typename R>
    sparse_polynomial<R> Z=sparse_polynomial<R>({0,1});

    template<typename R>
    sparse_polynomial<R> operator*(R a,const sparse_polynomial<R> &p)
    {
        auto q=p;
        return q*=a;
    }

}

#endif //CPLIBRARY_SPARSE_POLYNOMIAL_H

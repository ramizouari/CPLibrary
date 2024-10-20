//
// Created by ramizouari on 28/09/24.
//

#ifndef CPLIBRARY_LI_CHAO_H
#define CPLIBRARY_LI_CHAO_H
#include <vector>
#include <concepts>
#include <memory>
#include "algebra/abstract_algebra.h"
#include <bitset>

namespace cp::data_structures
{
    template<std::integral I, typename T>
    struct function : std::enable_shared_from_this<function<I,T>>
    {
        virtual T calculate(I n) = 0;
        T operator()(I n)
        {
            return calculate(n);
        }
    };


    template<std::integral I, typename T>
    struct translated_function : public function<I,T>
    {
        std::shared_ptr<function<I,T>> fn;
        T calculate(I n) override
        {
            return fn->calculate(n);
        }
    };



    template<std::integral I, typename T>
    struct li_chao_tree
    {

        using fn_ptr = std::shared_ptr<function<I,T>>;
        std::vector<std::vector<fn_ptr>> S;
        integer l;
        natural n,h;
        li_chao_tree(I l, I r) : l(l),n(std::bit_ceil<natural>(r-l))
        {
            int m=n;
            h=0;
            while(m)
            {
                m/=2;
                h++;
            }
            S.resize(h);
            for(int i=0;i<h;i++)
                S[i].resize(1<<i);

        }

        void add_function(fn_ptr fn)
        {

        }

        T minimum(I x)
        {
            return ;
        }

        bool compare(fn_ptr f, fn_ptr g, I x)
        {
            if(!g)
                return false;
            else if(!f)
                return true;
            else
                return f->calculate(x) < g->calculate(x);
        }

    private:

        void add_function(I l,I r,integer depth,fn_ptr f, std::bitset<2> s)
        {
            auto m=(l+r)/2;
            auto& g = S[depth][l >> (h-1-depth)];
            auto w= compare(f,g,m);
            if(w) g=f;
            if(s[0] ^ w)
            {
                auto s_ = s;
                s_[1]=r;
                add_function(l,m,depth+1,f,s_);
            }

            if(s[1] ^ r)
            {
                auto s_ = s;
                s_[0]=r;
                add_function(m,r,depth+1,f,s_);
            }

        }
    };
}
#endif //CPLIBRARY_LI_CHAO_H

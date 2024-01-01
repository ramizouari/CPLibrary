//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_DYN_FENWICK_TREE_H
#define CPLIBRARY_DYN_FENWICK_TREE_H
#include "algebra/binary_operation.h"
#include <memory>
namespace cp::data_structures::dynamic
{
    template<typename T>
    struct fenwick_tree {
        int n;
        invertible_binary_operation_ptr<T> F;
        std::vector<T> bit;

        fenwick_tree(int _n,std::shared_ptr<binary_operation<T>> _F, std::shared_ptr<invertible_operation<T>> _I):n(_n),F(_F,_I),bit(n,F.neutral_element()){}
        fenwick_tree(const std::vector<T> &X,std::shared_ptr<binary_operation<T>> _F, std::shared_ptr<invertible_operation<T>> _I) : fenwick_tree(X.size(),_F,_I)
        {
            for(int i=0;i<n;i++)
                update(i,X[i]);
        }
        T sum(int x) {
            if(x<0)
                return F.neutral_element();
            T ret = F.neutral_element();
            for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
                ret = F(ret,bit[i]);
            return ret;
        }

        T query(int a,int b)
        {
            return F(F.inv(sum(a-1)),sum(b));
        }

        T sum(int a,int b)
        {
            return query(a,b);
        }

        void add(int x, T delta) {
            for (int i = x; i < n; i = i | (i + 1))
                bit[i] = F(bit[i], delta);
        }

        void update(int x, T delta) {
            add(x,F(F.inv(sum(x,x)),delta));
        }
    };
}

#endif //CPLIBRARY_FENWICK_TREE_H

//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_FIXED_FENWICK_TREE_H
#define CPLIBRARY_FIXED_FENWICK_TREE_H
#include <vector>
namespace data_structures::fixed
{
    template<typename T,typename O>
    struct fenwick_tree {
        int n;
        std::vector<T> bit;
        inline static O F = O();

        fenwick_tree(int _n):n(_n),bit(n,O::neutral){}
        fenwick_tree(const std::vector<T> &X) : fenwick_tree(X.size())
        {
            for(int i=0;i<n;i++)
                update(i,X[i]);
        }
        T sum(int x) {
            if(x<0)
                return O::neutral;
            T ret = O::neutral;
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

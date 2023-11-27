//
// Created by ramizouari on 26/10/23.
//

#ifndef CPLIBRARY_FIXED_OTHERS_H
#define CPLIBRARY_FIXED_OTHERS_H
#include "data_structures/statistic_tree.h"
#include "segment_tree.h"
namespace cp::data_structures::fixed
{
    template<typename O>
    struct segment_matrix
    {
        using R=typename O::type;
        using type=R;
        std::vector<std::vector<segment_tree<O>>> S;
        std::vector<segment_tree<O>> segment_forest;
        std::vector<std::vector<R>> A;
        int n,h;
        segment_matrix(std::vector<std::vector<R>> &&_A):A(std::move(_A)),h(0)
        {
            n=bit_ceil(A.size());
            A.resize(n,std::vector<R>(bit_ceil(A[0].size()),O::neutral));
            int w=n;
            while(w)
            {
                w/=2;
                h++;
            }
            S.resize(h);
            int m=A[0].size();
            for(int i=0;i<m;i++)
            {
                std::vector<R> C;
                for(int j=0;j<n;j++)
                    C.push_back(A[j][i]);
                segment_forest.emplace_back(C);
            }
            std::vector<R> C(m);
            for(int i=h-1;i>=0;i--) for(int j=0;j<(1<<i);j++)
                {
                    for(int p=0;p<m;p++)
                        C[p]=segment_forest[p].S[i][j];
                    S[i].emplace_back(C);
                }
        }

        void update(int i,int j,R u)
        {
            A[i][j]=u;
            segment_forest[j].update(i,u);
            int r=h-1;
            while(r>=0)
            {
                S[r][i].update(j,segment_forest[j].S[r][i]);
                r--;
                i/=2;
            }
        }

        R query(int l,int r,int p,int q)
        {
            return query(l,r,p,q,0,n,0);
        }
    private:
        inline static O F=O();

        R query(int l,int r,int p,int q,int a,int b,int depth)
        {
            if(l>=r)
                return O::neutral;
            if(l==a && r==b)
                return S[depth][l>>(h-1-depth)].query(p,q);
            int mid=(a+b)/2;
            if(mid>r)
                return query(l,r,p,q,a,mid,depth+1);
            else if(mid<l)
                return query(l,r,p,q,mid,b,depth+1);
            else
                return F(query(l,mid,p,q,a,mid,depth+1),query(mid,r,p,q,mid,b,depth+1));
        }

    };

    template<typename T,typename O>
    struct fenwick_matrix {
        inline static O F=O();
        int n, m;
        std::vector<std::vector<T>> bit;

        fenwick_matrix(int _n,int _m):n(_n),m(_m),bit(n,std::vector<T>(m,O::neutral)){}
        int sum(int x, int y) {
            if(x<0||y<0)
                return O::neutral;
            int ret = 0;
            for (int i = x; i >= 0; i = (i & (i + 1)) - 1)
                for (int j = y; j >= 0; j = (j & (j + 1)) - 1)
                    ret = F(ret,bit[i][j]);
            return ret;
        }

        int sum(int a,int b,int c,int d)
        {
            //To Do
            //........................................
        }


        int query(int a,int b,int c,int d)
        {
            return F(F.inv(sum(a,c),sum(b,d)));
        }

        void add(int x, int y, int delta) {
            for (int i = x; i < n; i = i | (i + 1))
                for (int j = y; j < m; j = j | (j + 1))
                    bit[i][j] = F(bit[i][j],delta);
        }

        void update(int x, int y, int delta) {
            add(x,y,F(F.inv(sum(x,x,y,y)),delta));
        }
    };


    template<typename T,typename O>
    struct sparse_segment_tree
    {
        statistics_trees::sum_node<int, T, O>* tree;
    public:
        sparse_segment_tree():tree(nullptr) {}
        ~sparse_segment_tree()
        {
            destroy(tree);
            tree = nullptr;
        }

        void insert(int k, const T& v)
        {
            tree = statistics_trees::insert(tree, k, v);
        }

        void update(int k, const T& v)
        {
            tree = statistics_trees::insert_or_assign(k, v);
        }

        void erase(int k)
        {
            tree = statistics_trees::erase(tree, k);
        }

        T query(int l, int r) const
        {
            return sum(tree, l, r);
        }

        T index_query(int l, int r) const
        {
            return index_sum(tree,l, r);
        }
    };

    template<typename T, typename O>
    class dynamic_segment_tree
    {
        statistics_trees::sum_node<int, T, O>* tree;
        int size;
    public:
        dynamic_segment_tree() :tree(nullptr),size(0) {}
        ~dynamic_segment_tree()
        {
            destroy(tree);
            tree = nullptr;
            size = 0;
        }

        void push_back(const T& v)
        {
            tree = insert(tree, size++, v);
        }

        void update(int k, const T& v)
        {
            tree = insert_or_assign(k, v);
        }

        void pop_back()
        {
            tree = erase(tree, size--);
        }

        T query(int l, int r) const
        {
            return sum(tree, l, r);
        }

        T index_query(int l, int r) const
        {
            return query(l, r);
        }
    };

    template<typename T, typename O>
    struct ordered_segment_tree
    {
        statistics_trees::key_sum_node<T,O>* tree;
    public:
        ordered_segment_tree() :tree(nullptr) {}
        ~ordered_segment_tree()
        {
            destroy(tree);
            tree = nullptr;
        }

        void insert(const T& v)
        {
            tree = statistics_trees::insert(tree, v);
        }

        void erase(const T&v)
        {
            tree = statistics_trees::erase(tree, v);
        }

        T query(const T& l, const T& r) const
        {
            return key_sum(tree, l, r);
        }

        T index_query(int a, int b) const
        {
            return index_key_sum(tree, a, b);
        }
    };
}

#endif //CPLIBRARY_OTHERS_H

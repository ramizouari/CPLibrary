//
// Created by ramizouari on 17/11/2021.
//

#ifndef __STATISTIC_NODE__
#define __STATISTIC_NODE__
#include <tuple>
#include "binary_operation.h"


template<typename T,typename V,typename S>
struct statistic_node;

template<typename T,typename V,typename S>
struct statistic_node
{
    T v;
    V data;
    statistic_node *left,*right,*parent;
    int h;
    S statistic;
    statistic_node(T _v,V _data,statistic_node* _parent=nullptr):v(_v),data(_data),left(nullptr),right(nullptr),parent(_parent),h(1),statistic(v,data){}
    void update()
    {
        h=std::max(left?left->h:0,right?right->h:0)+1;
        S::update(this);
    }
};

template<typename T,typename V,typename S>
int height(statistic_node<T,V,S>*node)
{
    return node?node->h:0;
}

template<typename T,typename V,typename S>
int balance(statistic_node<T,V,S>* tree)
{
    return height(tree->left)-height(tree->right);
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* upper_bound(statistic_node<T,V,S>* tree,T v)
{
    if(tree==nullptr)
        return nullptr;
    if(tree->v<=v)
        return upper_bound(tree->right,v);
    else
    {
        auto w=upper_bound(tree->left,v);
        if(w==nullptr)
            return tree;
        return w;
    }
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* reverse_upper_bound(statistic_node<T,V,S>* tree,T v)
{
    if(tree==nullptr)
        return nullptr;
    if(tree->v>=v)
        return reverse_upper_bound(tree->left,v);
    else
    {
        auto w=reverse_upper_bound(tree->right,v);
        if(w==nullptr)
            return tree;
        return w;
    }
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* lower_bound(statistic_node<T,V,S>* tree,T v)
{
    if(tree==nullptr)
        return nullptr;
    if(tree->v<v)
        return lower_bound(tree->right,v);
    else
    {
        auto w=lower_bound(tree->left,v);
        if(w==nullptr)
            return tree;
        return w;
    }
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* reverse_lower_bound(statistic_node<T,V,S>* tree,T v)
{
    if(tree==nullptr)
        return nullptr;
    if(tree->v>v)
        return reverse_lower_bound(tree->left,v);
    else
    {
        auto w=reverse_lower_bound(tree->right,v);
        if(w==nullptr)
            return tree;
        return w;
    }
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* rebalance_right(statistic_node<T,V,S>* x)
{
    auto y=x->left,B=y->right;
    y->right=x;
    y->parent=x->parent;
    if(x->parent)
    {
        if(x->parent->left==x)
            x->parent->left=y;
        else x->parent->right=y;
    }
    x->parent=y;
    x->left=B;
    if(B) B->parent=x;
    x->update();
    return y;
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* rebalance_left(statistic_node<T,V,S>* x)
{
    auto y=x->right,B=y->left;
    y->left=x;
    y->parent=x->parent;
    if(x->parent)
    {
        if (x->parent->left == x)
            x->parent->left = y;
        else x->parent->right = y;
    }
    x->parent=y;
    x->right=B;
    if(B) B->parent=x;
    x->update();
    return y;
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* rebalance(statistic_node<T,V,S>* x)
{
    if(!x)
        return nullptr;
    if(balance(x)<-1)
    {
        if(balance(x->right)==1)
            rebalance_right(x->right);
        x=rebalance_left(x);
    }
    else if(balance(x)>1)
    {
        if(balance(x->left)==-1)
            rebalance_left(x->left);
        x=rebalance_right(x);
    }
    x->update();
    if(!x->parent)
        return x;
    return rebalance(x->parent);
}

template<typename T,typename V,typename S>
statistic_node<T,V,S>* insert(statistic_node<T,V,S>* tree,T v,V data,bool or_assign=false)
{
    if(!tree)
    {
        tree = new statistic_node<T,V,S>(v,data);
        return tree;
    }
    auto p=lower_bound(tree,v);
    if(p==nullptr)
    {
        p=tree;
        while(p->right) p=p->right;
    }
    else if (or_assign && p->v == v)
    {
        p->data = data;
        p->update();
        return rebalance(p);
    }
    else if(p->left)
    {
        p=p->left;
        while(p->right)
            p=p->right;
    }
    auto u=new statistic_node<T,V,S>(v,data,p);
    if(v<=p->v)
        p->left=u;
    else p->right=u;
    p->update();
    return rebalance(p);
}

template<typename T, typename V, typename S>
statistic_node<T, V, S>* insert_or_assign(statistic_node<T, V, S>* tree, T v, V data)
{
    return insert(tree, v, data, true);
}

template<typename T,typename V,typename S>
std::pair<statistic_node<T,V,S>*,statistic_node<T,V,S>*> extract(statistic_node<T,V,S>* tree,T v)
{
    auto p=lower_bound(tree,v);
    if(!p)
        return {nullptr,tree};
    if(!p->left)
    {
        auto w=p->parent?p->parent:p->right;
        if(p->parent)
        {
            if(p->parent->left==p) p->parent->left=p->right;
            else p->parent->right=p->right;
        }
        if(p->right) p->right->parent=p->parent;
        p->right=nullptr;
        p->parent=nullptr;
        return {p,rebalance(w)};
    }
    else if(!p->left->right)
    {
        auto w=p->left;
        if(p->parent)
        {
            if(p->parent->left==p) p->parent->left=p->left;
            else p->parent->right=p->left;
        }
        if(p->right) p->right->parent=w;
        w->right=p->right;
        w->parent=p->parent;
        p->right=nullptr;
        p->left=nullptr;
        p->parent=nullptr;
        return {p,rebalance(w)};
    }
    else
    {
        auto u=p->left;//Position of replacement
        while(u->right)
            u=u->right;
        auto s=u->parent;//Position of path to be updated
        s->right=u->left;
        if(u->left) u->left->parent=s;
        std::swap(u->v,p->v);
	    std::swap(u->data,p->data);
        u->left=nullptr;
        u->right=nullptr;
        u->parent=nullptr;
        return {u,rebalance(s)};
    }

}


template<typename T,typename V,typename S>
statistic_node<T,V,S>* erase(statistic_node<T,V,S>* tree,T v)
{
    auto P=extract(tree,v);
    delete P.first;
    return P.second;
}

template<typename T,typename S>
statistic_node<T,std::tuple<>,S>* insert(statistic_node<T,std::tuple<>,S>* tree,T v)
{
    return insert(tree,v,std::make_tuple());
}


template<typename T,typename V, typename S>
statistic_node<T,V,S>* update(statistic_node<T,V,S>*tree, T v, V data)
{
    auto p=lower_bound(tree,v);
    p->data=data;
    return rebalance(p);
}

template<typename T,typename V,typename S>
void destroy(statistic_node<T,V,S>*node)
{
    if(!node)
        return;
    destroy(node->left);
    destroy(node->right);
    delete node;
}

/*
* Order Statistic:
* It is a Self Balanced Binary Search Tree augmented with an order:
* - The order of an element can be calculated in O(log(n))
* - An element can be selected given its order in O(log(n))
*/

struct order_stats
{
    int size;
    order_stats() {}
    template<typename T,typename V>
    order_stats(T v,V data):size(1){}
    template<typename T,typename V>
    static void update(statistic_node<T,V,order_stats>*node);
};

template<typename T, typename V>
void order_stats::update(statistic_node<T, V, order_stats> *node) {
    node->statistic.size=(node->left?node->left->statistic.size:0)+1+(node->right?node->right->statistic.size:0);
}

template<typename T, typename V,typename OrderStats>
int tree_size(statistic_node<T, V, OrderStats> *node)
{
    return node?node->statistic.size:0;
}

template<typename T,typename V,typename OrderStats>
T order_inf(statistic_node<T,V, OrderStats> *tree,T v)
{
    if(!tree)
        return -1;
    if(v<tree->v)
        return order_inf(tree->left,v);
    else if(tree->v==v)
    {
        auto o=order_inf(tree->left,v);
        if(o!=-1)
            return o;
        else return tree_size(tree->left);
    }
    else
    {
        auto o=order_inf(tree->right,v);
        if(o!=-1)
            return tree_size(tree->left)+1+o;
        else return -1;
    }
}

template<typename T,typename V,typename OrderStats>
T order_sup(statistic_node<T,V, OrderStats> *tree,T v)
{
    if(!tree)
        return 0;
    if(v<tree->v)
        return order_sup(tree->left,v);
    else if(tree->v==v)
    {
        if(tree->right && tree->right->v==v)
            return tree_size(tree->left)+1+order_sup(tree->right,v);
        else return tree_size(tree->left);
    }
    else return tree_size(tree->left)+1+order_sup(tree->right,v);
}

template<typename T, typename V, typename OrderStats>
int order(statistic_node<T, V, OrderStats>* tree, T v)
{
    if (!tree)
        return 0;
    if (v < tree->v)
        return order(tree->left, v);
    else if (tree->v == v)
    {
        if (tree->right && tree->right->v == v)
            return tree_size(tree->left) + 1 + order(tree->right, v);
        else return tree_size(tree->left);
    }
    else return tree_size(tree->left) + 1 + order(tree->right, v);
}

template<typename T,typename V,typename OrderStats>
T select(statistic_node<T,V, OrderStats> *tree,int o)
{
    int s= tree_size(tree->left);
    if(s==o)
        return tree->v;
    else if(s<o)
        return select(tree->right,o-s-1);
    else return select(tree->left,o);
}

/*
* Sum Statistic:
* It is an Ordered Statistic Tree augmented with a sum acting on data:
* The sum is defined over an associative binary operation having a neutral element
* It supports range sum (L,R) for keys belonging to an interval [L,R[ 
*/
template<typename V, typename O>
struct sum_stats
{
    inline static constexpr O F = O();
    int size;
    V sum;
    sum_stats() {}
    template<typename T>
    sum_stats(T v, V data) :size(1), sum(data) {}
    template<typename T>
    static void update(statistic_node<T, V, sum_stats>* node);
    inline static const V& neutral = O::neutral;
};

template<typename V, typename O>
template<typename T>
void sum_stats<V, O>::update(statistic_node<T, V, sum_stats<V, O>>* node) {
    node->statistic.size = (node->left ? node->left->statistic.size : 0) + 1 + (node->right ? node->right->statistic.size : 0);
    node->statistic.sum = F(tree_sum(node->left),node->data, tree_sum(node->right));
}

template<typename T, typename V, typename SumStats>
V tree_sum(statistic_node<T, V, SumStats>* node)
{
    return node ? node->statistic.sum : SumStats::neutral;
}


template<typename T, typename V, typename SumStats>
V prefix_sum(statistic_node<T, V, SumStats>* tree, T U)
{
    if (!tree)
        return SumStats::neutral;
    else if (tree->v >= U)
        return prefix_sum(tree->left, U);
    else return SumStats::F(tree_sum(tree->left), tree->data, prefix_sum(tree->right, U));
}

template<typename T, typename V, typename SumStats>
V suffix_sum(statistic_node<T, V, SumStats>* tree, T L)
{
    if (!tree)
        return SumStats::neutral;
    else if (tree->v < L)
        return suffix_sum(tree->right, L);
    else return SumStats::F(suffix_sum(tree->left, L), tree->data, tree_sum(tree->right));
}

template<typename T, typename V, typename SumStats>
V sum(statistic_node<T, V, SumStats>* tree, T L, T R)
{
    if (!tree)
        return SumStats::neutral;
    if (tree->v < L)
        return sum(tree->right, L, R);
    else if (tree->v >= R)
        return sum(tree->left, L, R);
    else return SumStats::F(suffix_sum(tree->left, L), tree->data, prefix_sum(tree->right, R));
}

/*
* Key-sum Statistic
*/

template<typename T, typename O>
struct key_sum_stats
{
    inline static constexpr O F = O();
    int size;
    T key_sum;
    key_sum_stats() {}
    template<typename V>
    key_sum_stats(T v, V data) :size(1), key_sum(v) {}
    template<typename V>
    static void update(statistic_node<T, V, key_sum_stats>* node);
    inline static const T& key_neutral = O::neutral;
};

template<typename T, typename O>
template<typename V>
void key_sum_stats<T, O>::update(statistic_node<T, V, key_sum_stats<T, O>>* node) {
    node->statistic.size = (node->left ? node->left->statistic.size : 0) + 1 + (node->right ? node->right->statistic.size : 0);
    node->statistic.key_sum = F(tree_key_sum(node->left), node->v, tree_key_sum(node->right));
}


template<typename T, typename V, typename KeySumStats>
T tree_key_sum(statistic_node<T, V, KeySumStats>* node)
{
    return node ? node->statistic.key_sum : KeySumStats::key_neutral;
}

template<typename T, typename V, typename KeySumStats>
V prefix_key_sum(statistic_node<T, V, KeySumStats>* tree, T U)
{
    if (!tree)
        return KeySumStats::key_neutral;
    else if (tree->v >= U)
        return prefix_key_sum(tree->left, U);
    else return KeySumStats::F(tree_key_sum(tree->left), tree->v, prefix_key_sum(tree->right, U));
}

template<typename T, typename V, typename KeySumStats>
V suffix_key_sum(statistic_node<T, V, KeySumStats>* tree, T L)
{
    if (!tree)
        return KeySumStats::key_neutral;
    else if (tree->v < L)
        return suffix_key_sum(tree->right, L);
    else return KeySumStats::F(suffix_key_sum(tree->left, L), tree->v, tree_key_sum(tree->right));
}

template<typename T, typename V, typename KeySumStats>
T key_sum(statistic_node<T, V, KeySumStats>* tree, T L, T R)
{
    if (!tree)
        return KeySumStats::key_neutral;
    if (tree->v < L)
        return key_sum(tree->right, L, R);
    else if (tree->v >= R)
        return key_sum(tree->left, L, R);
    else return KeySumStats::F(suffix_key_sum(tree->left, L), tree->v, prefix_key_sum(tree->right, R));
}

#endif
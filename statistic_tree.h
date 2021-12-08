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
statistic_node<T,V,S>* insert(statistic_node<T,V,S>* tree,T v,V data)
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

struct order_stats
{
    int size;
    order_stats(){}
    template<typename T,typename V>
    order_stats(T v,V data):size(1){}
    template<typename T,typename V>
    static void update(statistic_node<T,V,order_stats>*node);
    template<typename T,typename V>
    static int tree_size(statistic_node<T,V,order_stats>*node);
};

template<typename T, typename V>
void order_stats::update(statistic_node<T, V, order_stats> *node) {
    node->statistic.size=(node->left?node->left->statistic.size:0)+1+(node->right?node->right->statistic.size:0);
}

template<typename T, typename V>
int order_stats::tree_size(statistic_node<T, V, order_stats> *node)
{
    return node?node->statistic.size:0;
}

template<typename T,typename V>
T order_inf(statistic_node<T,V,order_stats> *tree,T v)
{
    using statistic_type=order_stats;
    if(!tree)
        return -1;
    if(v<tree->v)
        return order_inf(tree->left,v);
    else if(tree->v==v)
    {
        auto o=order_inf(tree->left,v);
        if(o!=-1)
            return o;
        else return order_stats::tree_size(tree->left);
    }
    else
    {
        auto o=order_inf(tree->right,v);
        if(o!=-1)
            return statistic_type::tree_size(tree->left)+1+o;
        else return -1;
    }
}

template<typename T,typename V>
T order_sup(statistic_node<T,V,order_stats> *tree,T v)
{
    using statistic_type=order_stats;
    if(!tree)
        return 0;
    if(v<tree->v)
        return order_sup(tree->left,v);
    else if(tree->v==v)
    {
        if(tree->right && tree->right->v==v)
            return statistic_type::tree_size(tree->left)+1+order_sup(tree->right,v);
        else return statistic_type::tree_size(tree->left);
    }
    else return statistic_type::tree_size(tree->left)+1+order_sup(tree->right,v);
}

template<typename T,typename V>
T select(statistic_node<T,V,order_stats> *tree,int o)
{
    using statistic_type=order_stats;
    int s=order_stats::tree_size(tree->left);
    if(s==o)
        return tree->v;
    else if(s<o)
        return select(tree->right,o-s-1);
    else return select(tree->left,o);
}

#endif
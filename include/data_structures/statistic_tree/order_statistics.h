//
// Created by ramizouari on 19/12/23.
//

#ifndef CPLIBRARY_FIXED_ORDER_STATISTICS_H
#define CPLIBRARY_FIXED_ORDER_STATISTICS_H
#include "stats_tree_base.h"
namespace cp::data_structures::stats_trees {

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
        order_stats(const T& v,const V& data):size(1){}
        template<typename T,typename V>
        static void update(statistic_node<T,V,order_stats>*node);
    };

    template<typename T, typename V>
    void order_stats::update(statistic_node<T, V, order_stats> *node) {
        node->statistic.size=(node->left?node->left->statistic.size:0)+1+(node->right?node->right->statistic.size:0);
    }

    template<typename T, typename V,typename OrderStats>
    int size(statistic_node<T, V, OrderStats> *node)
    {
        return node?node->statistic.size:0;
    }

    template<typename T,typename V,typename OrderStats>
    T order_inf(statistic_node<T,V, OrderStats> *tree,const typename std::common_type<T>::type& v)
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
            else return size(tree->left);
        }
        else
        {
            auto o=order_inf(tree->right,v);
            if(o!=-1)
                return size(tree->left)+1+o;
            else return -1;
        }
    }

    template<typename T,typename V,typename OrderStats>
    T order_sup(statistic_node<T,V, OrderStats> *tree,const typename std::common_type<T>::type& v)
    {
        if(!tree)
            return 0;
        if(v<tree->v)
            return order_sup(tree->left,v);
        else if(tree->v==v)
        {
            if(tree->right && tree->right->v==v)
                return size(tree->left)+1+order_sup(tree->right,v);
            else return size(tree->left);
        }
        else return size(tree->left)+1+order_sup(tree->right,v);
    }

    template<typename T, typename V, typename OrderStats>
    int order(statistic_node<T, V, OrderStats>* tree, const typename std::common_type<T>::type& v )
    {
        if (!tree)
            return 0;
        if (v < tree->v)
            return order(tree->left, v);
        else if (tree->v == v)
        {
            if (tree->right && tree->right->v == v)
                return size(tree->left) + 1 + order(tree->right, v);
            else return size(tree->left);
        }
        else return size(tree->left) + 1 + order(tree->right, v);
    }

    template<typename T,typename V,typename OrderStats>
    statistic_node<T,V,OrderStats>* select(statistic_node<T,V, OrderStats> *tree,int o)
    {
        int s= size(tree->left);
        if(s==o)
            return tree;
        else if(s<o)
            return select(tree->right,o-s-1);
        else return select(tree->left,o);
    }

    template<typename T, typename V, typename OrderStats>
    statistic_node<T, V, OrderStats>* median(statistic_node<T, V, OrderStats>* tree)
    {
        return select(tree, (size(tree)-1) / 2);
    }

/*
* The following aliases gives the possibility to a much shorter and compact code
* For example:
* 1. sum_node_t<int,double,multiplies_t> is a shorthand for this atrocity:
*   statistic_node<int,double,sum_stats<double,multiplies_t<double>>>
*   this is a statistic_node whose key is of type int and values are doubles, the considered statistic is
*   an aggregation stats whose operation is multiplication
* The first two considers the case when the binary operator's class is a template
* The last two considers the case when the binary operator's class is not a template
*/

//

    template<typename T, typename V = std::monostate>
    using order_node = statistic_node<T, V, order_stats>;
}


#endif //CPLIBRARY_ORDER_STATISTICS_H

//
// Created by ramizouari on 19/12/23.
//

#ifndef CPLIBRARY_FIXED_VALUE_STATS_H
#define CPLIBRARY_FIXED_VALUE_STATS_H
#include "order_statistics.h"
namespace cp::data_structures::stats_trees
{

/**
* Sum Statistic:
* It is an Ordered Statistic Tree augmented with a sum acting on data:
* The sum is defined over an associative binary operation having a neutral element
* It supports range sum (L,R) for keys belonging to an interval [L,R[
* @Requirements
* 1. a class V
* 2. an associative binary operation O: VxV->V
* 3. O has a neutral element named 'neutral'. and it is defined as a public static attribute.
* 4. O has an overload for the ternary case: VxVxV->V, and its definition is compatible with the binary case.
* @Notes
* Formally, (V,O) is a monoid
*/
    template<typename V, typename O>
    struct sum_stats
    {
        inline static O F = O();
        int size;
        V sum;
        sum_stats() {}
        template<typename T>
        sum_stats(T v, V data) :size(1), sum(data) {}
        template<typename T>
        static void update(statistic_node<T, V, sum_stats>* node);
        inline static const V& neutral = O::neutral;
    };

    template<typename V>
    struct sum_stats<V,void>
    {
        int size;
        V sum;
        binary_operation_ptr<V> F;
        sum_stats() {}
        template<typename T>
        sum_stats(T v, V data,std::shared_ptr<binary_operation<V>> op=nullptr) :size(1), sum(data),F(op) {}
        template<typename T>
        static void update(statistic_node<T, V, sum_stats>* node);
        V neutral_element() const
        {
            return F->neutral_element();
        }
    };

    enum class subtree_type : bool
    {
        left,right
    };

    struct empty_tree_exception :public std::logic_error
    {
        empty_tree_exception():std::logic_error("tree is empty"){}
    };

    template<typename V, typename O>
    template<typename T>
    void sum_stats<V, O>::update(statistic_node<T, V, sum_stats<V, O>>* node) {
        node->statistic.size = (node->left ? node->left->statistic.size : 0) + 1 + (node->right ? node->right->statistic.size : 0);
        if constexpr (std::is_same_v<O, void>)
        {
            node->statistic.sum = node->statistic.F(tree_sum_dynamic(node, subtree_type::left), node->data,
                                                    tree_sum_dynamic(node, subtree_type::right));
        }
        else
            node->statistic.sum = F(tree_sum(node->left),node->data, tree_sum(node->right));

    }

    template<typename V>
    template<typename T>
    void sum_stats<V, void>::update(statistic_node<T, V, sum_stats<V, void>>* node) {
        node->statistic.size = (node->left ? node->left->statistic.size : 0) + 1 + (node->right ? node->right->statistic.size : 0);
        node->statistic.sum = node->statistic.F(tree_sum_dynamic(node, subtree_type::left), node->data,
                                                tree_sum_dynamic(node, subtree_type::right));

    }


    template<typename T,typename V>
    V tree_sum_dynamic(statistic_node<T, V, sum_stats<V, void>>* node, bool direction)
    {
        if (!node)
            throw empty_tree_exception();
        auto subtree=direction?node->right:node->left;
        return subtree?subtree->statistic.sum:node->statistic.F.neutral_element();
    }

    template<typename T,typename V>
    V tree_sum_dynamic(statistic_node<T, V, sum_stats<V, void>>* node, subtree_type direction)
    {
        return tree_sum_dynamic(node,static_cast<bool>(direction));
    }


    template<typename T, typename V, static_value_statistics<T,V> SumStats>
    V tree_sum(statistic_node<T, V, SumStats>* node)
    {
        return node ? node->statistic.sum : SumStats::neutral;
    }

    template<typename T, typename V, static_value_statistics<T,V> SumStats>
    V prefix_sum(statistic_node<T, V, SumStats>* tree,const
    typename std::common_type<T>::type& U)
    {
        if (!tree)
            return SumStats::neutral;
        else if (tree->v >= U)
            return prefix_sum(tree->left, U);
        else return SumStats::F(tree_sum(tree->left), tree->data, prefix_sum(tree->right, U));
    }

    template<typename T, typename V, dynamic_value_statistics<T,V> SumStats>
    V prefix_sum(statistic_node<T, V, SumStats>* tree,const
    typename std::common_type<T>::type& U)
    {
        if (!tree)
            throw empty_tree_exception();
        else if (tree->v >= U)
            return tree->left?prefix_sum(tree->left, U):tree->statistic.F.neutral_element();
        else return tree->statistic.F(tree_sum_dynamic(tree,subtree_type::left), tree->data, tree->right?prefix_sum(tree->right, U):tree->statistic.F.neutral_element());
    }

    template<typename T, typename V, static_value_statistics<T,V> SumStats>
    V suffix_sum(statistic_node<T, V, SumStats>* tree, const
    typename std::common_type<T>::type& L)
    {
        if (!tree)
            return SumStats::neutral;
        else if (tree->v < L)
            return suffix_sum(tree->right, L);
        else return SumStats::F(suffix_sum(tree->left, L), tree->data, tree_sum(tree->right));
    }

    template<typename T, typename V, dynamic_value_statistics<T,V> SumStats>
    V suffix_sum(statistic_node<T, V, SumStats>* tree, const
    typename std::common_type<T>::type& L)
    {
        if (!tree)
            throw empty_tree_exception();
        else if (tree->v < L)
            return tree->right?suffix_sum(tree->right, L):tree->statistic.F.neutral_element();
        else return tree->statistic.F(tree->left?suffix_sum(tree->left, L):tree->statistic.F.neutral_element(), tree->data, tree_sum_dynamic(tree,subtree_type::right));
    }

    template<typename T, typename V, static_value_statistics<T,V> SumStats>
    V sum(statistic_node<T, V, SumStats>* tree, const
    typename std::common_type<T>::type& L, const typename std::common_type<T>::type& R)
    {
        if (!tree)
            return SumStats::neutral;
        if (tree->v < L)
            return sum(tree->right, L, R);
        else if (tree->v >= R)
            return sum(tree->left, L, R);
        else return SumStats::F(suffix_sum(tree->left, L), tree->data, prefix_sum(tree->right, R));
    }

    template<typename T, typename V, dynamic_value_statistics<T,V> SumStats>
    V sum(statistic_node<T, V, SumStats>* tree, const
    typename std::common_type<T>::type& L, const typename std::common_type<T>::type& R)
    {
        if (!tree)
            throw empty_tree_exception();
        if (tree->v < L)
            return tree->right?sum(tree->right, L, R):tree->statistic.F.neutral_element();
        else if (tree->v >= R)
            return tree->left?sum(tree->left, L, R):tree->statistic.F.neutral_element();
        else return tree->statistic.F(tree->left?suffix_sum(tree->left, L):tree->statistic.F.neutral_element(), tree->data,
                                      tree->right?prefix_sum(tree->right, R):tree->statistic.F.neutral_element());

    }

    template<typename T, typename V, static_value_statistics<T,V> SumStats>
    V prefix_index_sum(statistic_node<T, V, SumStats>* tree, int n)
    {
        if (!tree || n <= 0)
            return SumStats::neutral;
        else if (size(tree->left) >= n)
            return prefix_index_sum(tree->left, n);
        else return SumStats::F(tree_sum(tree->left), tree->data, prefix_index_sum(tree->right, n - 1 - size(tree->left)));

    }

    template<typename T, typename V, dynamic_value_statistics<T,V> SumStats>
    V prefix_index_sum(statistic_node<T, V, SumStats>* tree, int n)
    {
        if (!tree)
            throw empty_tree_exception();
        else if (size(tree->left) >= n)
            return tree->left?prefix_index_sum(tree->left, n):tree->statistic.F.neutral_element();
        else return tree->statistic.F(tree_sum_dynamic(tree,subtree_type::left), tree->data,
                                      tree->right?prefix_index_sum(tree->right, n - 1 - size(tree->left)) : tree->statistic.F.neutral_element());
    }

    template<typename T, typename V, static_value_statistics<T,V> SumStats>
    V suffix_index_sum(statistic_node<T, V, SumStats>* tree, int n)
    {
        if (!tree)
            return SumStats::neutral;
        else if (size(tree->right) >= n)
            return suffix_index_sum(tree->right, n);
        else return SumStats::F(suffix_index_sum(tree->left, n - 1 - size(tree->right)), tree->data, tree_sum(tree->right));
    }

    template<typename T, typename V, dynamic_value_statistics<T,V> SumStats>
    V suffix_index_sum(statistic_node<T, V, SumStats>* tree, int n)
    {
        if (!tree)
            throw empty_tree_exception();
        else if (size(tree->right) >= n)
            return tree->right?suffix_index_sum(tree->right, n):tree->statistic.F.neutral_element();
        else return tree->statistic.F(tree->left?suffix_index_sum(tree->left, n - 1 - size(tree->right)):tree->statistic.F.neutral_element(),
                                      tree->data,
                                      tree_sum_dynamic(tree->right));
    }

    template<typename T, typename V, typename SumStats>
    V index_sum(statistic_node<T, V, SumStats>* tree, int a, int b)
    {
        if (!tree || a >= b)
            return SumStats::neutral;
        int left_size = size(tree->left);
        if (left_size < a)
            return SumStats::F(tree->data,index_sum(tree->right, a - left_size - 1, b - left_size - 1));
        else if (left_size >= b)
            return index_sum(tree->left, a, b);
        else return SumStats::F(suffix_index_sum(tree->left, left_size), tree->data, prefix_index_sum(tree->right, left_size - 1));
    }

    template<typename T, typename V, template<typename S = T> typename O>
    using sum_node_t = statistic_node<T, V, sum_stats<V, O<V>>>;

    template<typename T, typename V, typename O>
    using sum_node = statistic_node<T, V, sum_stats<V, O>>;

    template<typename T, typename V>
    using dynamic_sum_node = statistic_node<T, V, sum_stats<V, void>>;

}

#endif //CPLIBRARY_VALUE_STATS_H

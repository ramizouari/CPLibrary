//
// Created by ramizouari on 19/12/23.
//

#ifndef CPLIBRARY_FIXED_KEY_STATISTICS_H
#define CPLIBRARY_FIXED_KEY_STATISTICS_H
#include "order_statistics.h"
namespace cp::data_structures::stats_trees
{

/** Key Sum Statistic:
* It is an Ordered Statistic Tree augmented with a sum acting on keys:
* The sum is defined over an associative binary operation having a neutral element
* It supports range sum (L,R) for keys belonging to the interval [L,R[
* @Requirements
* 1. a class T
* 2. an associative binary operation O: TxT->T
* 3. O has a neutral element named 'neutral'. and it is defined as a public static attribute.
* 4. O has an overload for the ternary case: TxTxT->T, and its definition is compatible with the binary case.
* @Notes
* Formally, (T,O) is a monoid
* The order of T does not need to be compatible with O
*/

    template<typename T, typename O>
    struct key_sum_stats
    {
        inline static O F = O();
        int size;
        T key_sum;
        key_sum_stats() {}
        template<typename V>
        key_sum_stats(T v, V data) :size(1), key_sum(v) {}
        template<typename V>
        static void update(statistic_node<T, V, key_sum_stats>* node);
        inline static const T& key_neutral = O::neutral;
    };

    template<typename T>
    struct key_sum_stats<T,void>
    {
        int size;
        T key_sum;
        binary_operation_ptr<T> F;
        key_sum_stats() {}
        template<typename V>
        key_sum_stats(T v, V data) :size(1), key_sum(v) {}
        template<typename V>
        static void update(statistic_node<T, V, key_sum_stats>* node);
        T key_neutral_element() const { return F->neutral_element(); }
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
    T prefix_key_sum(statistic_node<T, V, KeySumStats>* tree, const
    typename std::common_type<T>::type &U)
    {
        if (!tree)
            return KeySumStats::key_neutral;
        else if (tree->v >= U)
            return prefix_key_sum(tree->left, U);
        else return KeySumStats::F(tree_key_sum(tree->left), tree->v, prefix_key_sum(tree->right, U));
    }

    template<typename T, typename V, typename KeySumStats>
    T suffix_key_sum(statistic_node<T, V, KeySumStats>* tree, const
    typename std::common_type<T>::type &L)
    {
        if (!tree)
            return KeySumStats::key_neutral;
        else if (tree->v < L)
            return suffix_key_sum(tree->right, L);
        else return KeySumStats::F(suffix_key_sum(tree->left, L), tree->v, tree_key_sum(tree->right));
    }

    template<typename T, typename V, typename KeySumStats>
    T key_sum(statistic_node<T, V, KeySumStats>* tree, const typename std::common_type<T>::type& L,
              const typename std::common_type<T>::type& R)
    {
        if (!tree)
            return KeySumStats::key_neutral;
        if (tree->v < L)
            return key_sum(tree->right, L, R);
        else if (tree->v >= R)
            return key_sum(tree->left, L, R);
        else return KeySumStats::F(suffix_key_sum(tree->left, L), tree->v, prefix_key_sum(tree->right, R));
    }

    template<typename T, typename V, typename KeySumStats>
    T prefix_index_key_sum(statistic_node<T, V, KeySumStats>* tree, int n)
    {
        if (!tree || n <= 0)
            return KeySumStats::key_neutral;
        else if (size(tree->left) >= n)
            return prefix_index_key_sum(tree->left, n);
        else return KeySumStats::F(tree_key_sum(tree->left), tree->v, prefix_index_key_sum(tree->right, n - 1 - size(tree->left)));
    }

    template<typename T, typename V, typename KeySumStats>
    T suffix_index_key_sum(statistic_node<T, V, KeySumStats>* tree, int n)
    {
        if (!tree)
            return KeySumStats::key_neutral;
        else if (size(tree->right) >= n)
            return suffix_index_key_sum(tree->right, n);
        else return KeySumStats::F(suffix_index_key_sum(tree->left, n-1-size(tree->right)), tree->v, tree_key_sum(tree->right));
    }

    template<typename T, typename V, typename KeySumStats>
    T index_key_sum(statistic_node<T, V, KeySumStats>* tree, int a,int b)
    {
        if (!tree || a>=b)
            return KeySumStats::key_neutral;
        int left_size = size(tree->left);
        if (left_size < a)
            return KeySumStats::F(tree->v,index_key_sum(tree->right, a-left_size-1, b-left_size-1));
        else if (left_size >= b)
            return index_key_sum(tree->left, a, b);
        else return KeySumStats::F(suffix_index_key_sum(tree->left, left_size-a), tree->v, prefix_index_key_sum(tree->right, b-left_size-1));
    }

    template<typename T, template<typename S = T> typename O, typename V=std::monostate>
    using key_sum_node_t = statistic_node<T, V, key_sum_stats<T, O<T>>>;

    template<typename T, typename O, typename V = std::monostate>
    using key_sum_node = statistic_node<T, V, key_sum_stats<T, O>>;

}
#endif //CPLIBRARY_KEY_STATISTICS_H

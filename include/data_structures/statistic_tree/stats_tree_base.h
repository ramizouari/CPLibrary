//
// Created by ramizouari on 19/12/23.
//
#ifndef CPLIBRARY_STATISTICS_TREE_BASE_H
#define CPLIBRARY_STATISTICS_TREE_BASE_H
#include <tuple>
#include <memory>
#include "algebra/binary_operation.h"
#include <variant>
namespace cp::data_structures::stats_trees
{
    /**
* @definition Order Statistic Tree: It is an AVL-tree augmented with a statistic.
* The statistic is a function of:
* </ol>
* <li> The key
* <li> The value
* <li> left and right subtree
* </ol>
* @Requirements
* <ol><li> a class T
* <li> (T,<=) is a totally ordered set
* <li> S is a statistic type
* </ol>
* @tparam T the type of the key
* @tparam V the type of the value
* @tparam S the statistic type: a statistic S is a class having some attributes serving as additional informations
* Those informations can be aggregated via the static update method
* @Requirements
* <ol>
* <li> class S
* <li> S has a public static method called update which accepts Tree<T,V,S>.
* <li> the update method "updates" adequately the statistics, and only the statistics
* </ol>
*/


    template<typename T,typename V,typename S>
    struct statistic_node;

    template<typename Stats,typename T,typename V>
    concept base_statistics = requires(Stats s,T key, V data,statistic_node<T,V,Stats>* node)
    {
        {Stats::update(node)}->std::same_as<void>;
        Stats();
        Stats(key,data);
        {s.size} -> std::convertible_to<int>;
    };
    template<typename Stats,typename T,typename V>
    concept dynamic_value_statistics = base_statistics<Stats,T,V> && requires(Stats s,T key, V data)
    {
        {s.neutral_element()}->std::convertible_to<V>;
        {s.F} -> std::convertible_to<binary_operation_ptr<V>>;
    };
    template<typename Stats,typename T,typename V>
    concept dynamic_key_statistics = base_statistics<Stats,T,V> && requires(Stats s,T key, V data)
    {
        Stats::Stats(key,data);
        {s.key_neutral_element()}->std::convertible_to<T>;
        {s.F} -> std::convertible_to<binary_operation_ptr<T>>;
    };

    template<typename Stats,typename T,typename V>
    concept static_value_statistics = base_statistics<Stats,T,V> && requires(Stats s,T key, V data)
    {
        {Stats::F};
        {Stats::neutral} -> std::convertible_to<V>;
    };

    template<typename T,typename V,typename S>
    struct statistic_node
    {
        T v;
        V data;
        int h;
        S statistic;
        statistic_node *left,*right,*parent;
        statistic_node(T _v,V _data,statistic_node* _parent=nullptr):v(_v),data(_data),left(nullptr),right(nullptr),parent(_parent),h(1),statistic(v,data)
        {
            if constexpr (dynamic_key_statistics<S,T,V> || dynamic_value_statistics<S,T,V>)
            {
                if(parent)
                    statistic.F=parent->statistic.F;
            }
        }

        void update()
        {
            h=std::max(left?left->h:0,right?right->h:0)+1;
            S::update(this);
        }
    };

/*
* Get the height of the (sub)-tree
*/
    template<typename T,typename V,typename S>
    int height(statistic_node<T,V,S>*node)
    {
        return node?node->h:0;
    }

/*
* Get the balance of the current node
*/
    template<typename T,typename V,typename S>
    int balance(statistic_node<T, V, S>* tree)
    {
        return height(tree->left) - height(tree->right);
    }

/**
* Find node that has the given key, return nullptr otherwise
* @Notes
* 1. If the value is altered, and the statistic depends on that value. The calculated statistic will not be updated.
*   In that case, better use the insert_or_assign function or the update function
* 2. If the key is altered, It will probably cause the tree to be in an inconsistent state unless the tree does not have a key
* on the interval whose limit points are the old key value, and the new one.
*/
    template<typename T, typename V, typename S>
    statistic_node<T,V,S>* find(statistic_node<T, V, S>* node, const typename std::common_type<T>::type& v)
    {
        if (!node)
            return nullptr;
        if (node->v == v)
            return node->data;
        else if (node->v < v)
            return find(node->right, v);
        else return find(node->left, v);
    }

/*
* Find the data mapped by the given key
* @Requirements
* The Order Statistic Tree must contain at least one such key
* @Exception
* std::logic_error thrown if no such key is found
* @Notes
* If the value is altered, and the statistic depends on that value. The calculated statistic will not be updated.
* In that case, better use the insert_or_assign function or the update function
*/
    template<typename T, typename V, typename S>
    V& value_at(statistic_node<T, V, S>* node, const typename std::common_type<T>::type& v)
    {
        [[unlikely]]
        if (!node)
            throw std::out_of_range("key does not exist");
        if (node->v == v)
            return node;
        else if (node->v < v)
            return value_at(node->right, v);
        else return value_at(node->left, v);
    }

/*
* Get the node whose key is strictly bigger than v
*/
    template<typename T,typename V,typename S>
    statistic_node<T,V,S>* upper_bound(statistic_node<T,V,S>* tree,const
    typename std::common_type<T>::type& v)
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

/*
* Get the node whose key is strictly smaller than v
*/
    template<typename T,typename V,typename S>
    statistic_node<T,V,S>* reverse_upper_bound(statistic_node<T,V,S>* tree,
                                               const typename std::common_type<T>::type & v)
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

/*
* Get the node whose key is not smaller than v
*/
    template<typename T,typename V,typename S>
    statistic_node<T,V,S>* lower_bound(statistic_node<T,V,S>* tree,
                                       const typename std::common_type<T>::type& v)
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
    statistic_node<T,V,S>* reverse_lower_bound(statistic_node<T,V,S>* tree,
                                               const typename std::common_type<T>::type &v)
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

/*
* Get the node with the smallest key
*/
    template<typename T, typename V, typename S>
    statistic_node<T, V, S>* begin(statistic_node<T, V, S>* tree)
    {
        if (!tree)
            return nullptr;
        while (tree->left)
            tree = tree->left;
        return tree;
    }

/*
* Get the successor of the current node
*/
    template<typename T, typename V, typename S>
    statistic_node<T, V, S>* next(statistic_node<T, V, S>* tree)
    {
        if (tree == nullptr)
            return nullptr;
        if (tree->right)
        {
            tree = tree->right;
            while (tree->left)
                tree = tree->left;
            return tree;
        }
        else
        {
            auto tmp = tree;
            tree = tree->parent;
            while (tree && tree->v < tmp->v)
                tree = tree->parent;
            if (!tree)
                return nullptr;
            return tree;
        }
    }

/*
* Get the previous of the current node
*/
    template<typename T, typename V, typename S>
    statistic_node<T, V, S>* prev(statistic_node<T, V, S>* tree)
    {
        if (tree == nullptr)
            return nullptr;
        if (tree->left)
        {
            tree = tree->left;
            while (tree->right)
                tree = tree->right;
            return tree;
        }
        else
        {
            auto tmp = tree;
            tree = tree->parent;
            while (tree && tree->v > tmp->v)
                tree = tree->parent;
            if (!tree)
                return nullptr;
            return tree;
        }
    }

/*
* Applies a right rotation to the ordered statistic tree on the current node.
*/
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

/*
* Applies a left rotation to the ordered statistic tree on the current node.
*/
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

/*
* Rebalance the ordered statistic tree.
*/
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

/*
* Insert (v,data) into the ordered statistic tree
* @Cases
* 1. If or_assign is false, (v,data) is inserted while allowing duplicates
* 2. if or_assign is true, (v,data) is inserted if the tree does not have a key v. Otherwise the value mapped by v is changed to data.
*/
    template<typename T,typename V,typename S>
    statistic_node<T,V,S>* insert(statistic_node<T,V,S>* tree,const typename std::common_type<T>::type& v,
                                  const typename std::common_type<V>::type& data,bool or_assign=false)
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

    template<typename T, typename S>
    statistic_node<T, std::monostate, S>* insert(statistic_node<T, std::monostate, S>* tree,
                                                 const typename std::common_type<T>::type& v,bool or_assign=false)
    {
        return insert(tree, v, {},or_assign);
    }

    template<typename T,typename V,typename S,typename BinaryOpPtr>
    statistic_node<T,V,S>* insert(statistic_node<T,V,S>* tree,const typename std::common_type<T>::type& v,
                                  const typename std::common_type<V>::type& data, BinaryOpPtr F,bool or_assign=false)
    {
        tree=insert(tree,v,data,or_assign);
        tree->statistic.F=F;
        return tree;
    }

/*
* Insert (v,data) into the ordered statistic tree if it does not have a key v
* Otherwise, change the value mapped by v to data.
*/
    template<typename T, typename V, typename S>
    statistic_node<T, V, S>* insert_or_assign(statistic_node<T, V, S>* tree,const
    typename std::common_type<T>::type& v,const typename std::common_type<V>::type& data)
    {
        return insert(tree, v, data, true);
    }

    template<typename T, typename V, typename S,typename BinaryOpPtr>
    statistic_node<T, V, S>* insert_or_assign(statistic_node<T, V, S>* tree,const
    typename std::common_type<T>::type& v,const typename std::common_type<V>::type& data, BinaryOpPtr F)
    {
        return insert(tree, v, data,F, true);
    }

/*
* Insert (v,data) into the ordered statistic tree if it does not have a key v
* Otherwise, Do nothing
*/
    template<typename T, typename S>
    statistic_node<T, std::monostate, S>* insert_or_assign(statistic_node<T, std::monostate, S>* tree, const
    typename std::common_type<T>::type& v)
    {
        return insert_or_assign(tree, v, {});
    }

/*
* Extract a node from an ordered statistic tree given its key
* The node is not deleted
* @Requirements
* The tree does have a node with the given key
*/
    template<typename T,typename V,typename S>
    std::pair<statistic_node<T,V,S>*,statistic_node<T,V,S>*> extract(statistic_node<T,V,S>* tree,const
    typename std::common_type<T>::type& v)
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
    statistic_node<T,V,S>* erase(statistic_node<T,V,S>* tree,
                                 const typename std::common_type<T>::type& v)
    {
        auto P=extract(tree,v);
        delete P.first;
        return P.second;
    }

    template<typename T,typename V, typename S>
    statistic_node<T,V,S>* update(statistic_node<T,V,S>*tree,
                                  const typename std::common_type<T>::type&  v,
                                  const typename std::common_type<V>::type& data)
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
* Merge a node with two trees
* @Requirements
* 1. all keys of left are less or equal to the key of the root node
* 2. all keys of right are greater or equal to the key of the root node
* 3. the left & right trees do not have a parent
*/
    template<typename T,typename V,typename S>
    statistic_node<T, V, S>* merge_with_root(statistic_node<T, V, S>* root,
                                             statistic_node<T, V, S> *left, statistic_node<T, V, S> *right)
    {
        if (left == nullptr)
            return insert(right, root->v, root->data);
        else if (right == nullptr)
            return insert(left, root->v, root->data);
        auto potential = height(left) - height(right);
        if (potential <=0)
        {
            while (potential < -1)
                potential = height(left) - height(right=right->left);
            if (right && right->parent)
                right->parent->left = root;
            if(right)
                root->parent = right->parent;
        }
        else
        {
            while (potential > 1)
                potential = height(left=left->right) - height(right);
            if (left && left->parent)
                left->parent->right = root;
            if(left)
                root->parent = left->parent;
        }
        root->left = left;
        root->right = right;
        root->update();
        if (left)
            left->parent = root;
        if (right)
            right->parent = root;
        return rebalance(root);
    }

/*
* Merge two trees
* @Requirements
* The biggest key of left is smaller than the smallest key of right
*/
    template<typename T, typename V, typename S>
    statistic_node<T, V, S>* merge(statistic_node<T, V, S>* left, statistic_node<T, V, S>* right)
    {
        if (!left)
            return right;
        statistic_node<T, V, S>* last = left;
        while (last->right)
            last = last->right;
        auto [root,L] = extract(last, last->v);
        return merge_with_root(root, L, right);
    }
/*
* Split a tree into two trees with respect to threshold
* - The left part is the resultant tree whose keys are smaller than threshold
* - The right part is the resultant tree whose keys are not smaller than threshold
*/
    template<typename T, typename V, typename S>
    std::pair<statistic_node<T, V, S>*,statistic_node<T,V,S>*>
    split(statistic_node<T, V, S>* node,T threshold)
    {
        statistic_node<T, V, S>* left = nullptr, * right = nullptr;
        if (!node)
            return std::make_pair(left, right);
        if (node->right)
            node->right->parent = nullptr;
        if (node->left)
            node->left->parent = nullptr;
        if (node->v < threshold)
        {
            auto [L, R] = split(node->right, threshold);
            if (L)
                L->parent = nullptr;
            if (R)
                R->parent = nullptr;
            left = merge_with_root(node,node->left,L);
            right = R;
        }
        else
        {
            auto [L, R] = split(node->left, threshold);
            if (L)
                L->parent = nullptr;
            if (R)
                R->parent = nullptr;
            right = merge_with_root(node,R, node->right);
            left = L;
        }
        return std::make_pair(left, right);
    }
}
#endif //CPLIBRARY_ORDER_STATISTICS_H

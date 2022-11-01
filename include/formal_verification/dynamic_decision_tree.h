//
// Created by ramizouari on 18/10/22.
//

#ifndef VFSF_DYNAMIC_DECISION_TREE_H
#define VFSF_DYNAMIC_DECISION_TREE_H
#include "abstract_decision_tree.h"
#include <unordered_map>
#include <memory>
#include <functional>

namespace  FormalSpecification
{
    struct dynamic_decision_tree : public abstract_decision_tree
    {
        std::vector<int> permutation;
        int n;
        std::unordered_map<node_t,node_t*> representative;
        explicit dynamic_decision_tree(int n);
        explicit dynamic_decision_tree(std::vector<int> _permutation);
        dynamic_decision_tree(const dynamic_decision_tree &O);
        dynamic_decision_tree(dynamic_decision_tree&&O) noexcept = default;
        ~dynamic_decision_tree() override;
        virtual bool evaluate(const std::vector<bool> &S) const =0;

        bool evaluate_tree(const std::vector<bool> &S) const;

        bool operator()(const std::vector<bool> &S) const;

        template<typename ...T>
        bool operator()(T ...x)
        {
            std::vector<bool> S{x...};
            return evaluate(S);
        }

        void build();

        natural states_count() const;

        natural expanded_states_count() const;

        static void reseed(natural seed);


        dynamic_decision_tree& invert();


    protected:

        void compress_join();

        /**
         * @brief This function applies the first rule.
         * @details It transforms the complete tree to a directed acyclic graph by fusing equivalent states into one state.
         * Equivalent states are the states whose underlying subtrees are isomorphic
         * @param N The node on which the direct children will be processed
         *
         * */
        void compress(node_t &N);

        /**
         * @brief This function applies the second rule.
         * @details It eliminates recursively redundant states.
         * @note A redundant state is a state which has equivalent (same) successors
         * @warning This function invalidates the 'representative' map. After calling this function, we should recalculate it.
         * */
        void join();

        /**
        * @brief This function applies the second rule.
        * @details It eliminates recursively redundant states.
        * @note A redundant state is a state which has equivalent (same) successors
        * @warning This function invalidates the 'representative' map. After calling this function, we should recalculate it.
        * @param N the state who has potential redundant children
        * @param mapper a mapping between deleted redundant states and their replacement
        * */
        void join(node_t&N,std::unordered_map<node_t*,node_t*> &mapper);

        void recalculate_representative(node_t*N);

        node_t* create_left(std::vector<bool> &S,int depth);

        node_t* create_right(std::vector<bool> &S,int depth);

        node_t* build(std::vector<bool> &S,int depth=0);

        friend class default_dynamic_decision_tree operator~(const dynamic_decision_tree &A);

    };


    struct default_dynamic_decision_tree : public dynamic_decision_tree
    {
        using dynamic_decision_tree::dynamic_decision_tree;
        bool evaluate(const std::vector<bool> &S) const override;

        template<typename BooleanOperator>
        static default_dynamic_decision_tree merge(const dynamic_decision_tree &u,const dynamic_decision_tree &v , BooleanOperator && B)
        {
            if(u.permutation!=v.permutation)
                throw std::runtime_error("Permutations must be equal");
            std::vector<int> O(u.permutation.size());
            for(int i=0;i<O.size();i++)
                O[u.permutation[i]]=i;
            default_dynamic_decision_tree X(u.permutation);
            X.root=FormalSpecification::merge(u.root,v.root,O,std::forward<BooleanOperator>(B));
            X.compress_join();
            return X;
        }
    };

    struct table_based_decision_tree : public  dynamic_decision_tree
    {
        std::vector<bool> T;
        table_based_decision_tree(int n,const std::vector<bool>& _T);
        table_based_decision_tree(int n,std::vector<bool> && _T);
        table_based_decision_tree(const std::vector<int> &permutation,const std::vector<bool>& _T);
        table_based_decision_tree(std::vector<int> &&permutation,const std::vector<bool>& _T);
        table_based_decision_tree(const std::vector<int> &permutation, std::vector<bool>&& _T);
        table_based_decision_tree(std::vector<int> &&permutation, std::vector<bool>&& _T);
        bool evaluate(const std::vector<bool> &S) const override;
    };


    struct lambda_based_decision_tree : public  dynamic_decision_tree
    {
        std::function<bool(std::vector<bool>)> T;
        lambda_based_decision_tree(int n,const std::function<bool(std::vector<bool>)> & _T);
        lambda_based_decision_tree(int n, std::function<bool(std::vector<bool>)> && _T);
        lambda_based_decision_tree(const std::vector<int> &permutation,const std::function<bool(std::vector<bool>)>& _T);
        lambda_based_decision_tree(const std::vector<int> &permutation,std::function<bool(std::vector<bool>)>&& _T);
        lambda_based_decision_tree(std::vector<int> &&permutation,const std::function<bool(std::vector<bool>)>& _T);
        lambda_based_decision_tree(std::vector<int> &&permutation, std::function<bool(std::vector<bool>)>&& _T);

        bool evaluate(const std::vector<bool> &S) const override;
    };


    template<typename BooleanOperator>
    default_dynamic_decision_tree merge(const dynamic_decision_tree &u,const dynamic_decision_tree &v , BooleanOperator && B)
    {
        return default_dynamic_decision_tree::merge(u,v,std::forward<BooleanOperator>(B));
    }


    default_dynamic_decision_tree operator&(const dynamic_decision_tree& A,const dynamic_decision_tree &B);
    default_dynamic_decision_tree operator|(const dynamic_decision_tree& A,const dynamic_decision_tree &B);
    default_dynamic_decision_tree operator^(const dynamic_decision_tree& A,const dynamic_decision_tree &B);
    default_dynamic_decision_tree operator~(const dynamic_decision_tree &A);
}




#endif //VFSF_DYNAMIC_DECISION_TREE_H

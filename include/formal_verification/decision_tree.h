//
// Created by ramizouari on 18/10/22.
//

#ifndef VFSF_DECISION_TREE_H
#define VFSF_DECISION_TREE_H
#include <bitset>
#include <random>
#include <array>
#include <unordered_map>
#include "abstract_decision_tree.h"

namespace FormalSpecification
{
    template<int n>
    struct decision_tree : public abstract_decision_tree
    {
        node_t* root;
        std::array<int,n> permutation;
        std::unordered_map<node_t,node_t*> representative;
        explicit decision_tree(): abstract_decision_tree(n)
        {
            for(int i=0;i<n;i++)
                permutation[i]=i;
        }
        explicit decision_tree(std::array<int,n> permutation):permutation(std::move(permutation)),
                                                              abstract_decision_tree(n)
        {
        }

        decision_tree(decision_tree &O): abstract_decision_tree(n)
        {
            root=node_t::clone(O.root);
            permutation=O.permutation;
            recalculate_representative(root);
        }

        decision_tree(decision_tree &&O)  noexcept =default;

        virtual ~decision_tree()
        {
            for(auto &[_,node]:representative)
                delete node;
        };

        virtual bool evaluate(const std::bitset<n> &S) const =0;

        bool evaluate_tree(const std::bitset<n> &S) const
        {
            auto node=root;
            while(node->symbol<n)
            {
                if(S[node->symbol])
                    node=node->right;
                else
                    node=node->left;
            }
            return static_cast<Truth>(node->symbol)==Truth::True;
        }

        bool operator()(const std::bitset<n> &S) const
        {
            return evaluate(S);
        }

        template<typename ...T>
        bool operator()(T ...x) requires (sizeof...(T) == n)
        {
            std::bitset<n> S{x...};
            return evaluate(S);
        }

        void build()
        {
            std::bitset<n> S;
            root=build(S);
            compress_join();
        }

        natural states_count() const
        {
            return representative.size();
        }

        natural expanded_states_count() const
        {
            return (1<<(n+1))-1;
        }


        static void reseed(natural seed)
        {
            node_t::g.seed(seed);
        }

    private:

        void compress_join()
        {
            representative[*root]=root;
            compress(*root);
            join();
            representative.clear();
            recalculate_representative(root);
        }

        /**
         * @brief This function applies the first rule.
         * @details It transforms the complete tree to a directed acyclic graph by fusing equivalent states into one state.
         * Equivalent states are the states whose underlying subtrees are isomorphic
         * @param N The node on which the direct children will be processed
         *
         * */
        void compress(node_t &N)
        {
            if(N.left)
            {
                //Apply a recursive call
                compress(*N.left);
                //If The left child is isomorphic to some previously seen tree, replace it with the representative
                if(representative.contains(*N.left) && representative[*N.left]!=N.left)
                {
                    auto tmp=N.left;
                    N.left=representative[*N.left];
                    delete tmp;
                }
                else
                    representative[*N.left]=N.left;
            }

            if(N.right)
            {
                compress(*N.right);
                if (representative.contains(*N.right) && representative[*N.right] != N.right) {
                    auto tmp = N.right;
                    N.right = representative[*N.right];
                    delete tmp;
                } else
                    representative[*N.right] = N.right;

            }
        }

        /**
     * @brief This function applies the second rule.
     * @details It eliminates recursively redundant states.
     * @note A redundant state is a state which has equivalent (same) successors
     * @warning This function invalidates the 'representative' map. After calling this function, we should recalculate it.
     * */
        void join()
        {
            std::unordered_map<node_t*,node_t*> mapper;
            join(*root,mapper);
            //If the root is redundant, remove it
            if(root->left && root->left == root->right)
            {
                auto tmp=root;
                root=root->left;
                delete tmp;
            }
        }

        /**
    * @brief This function applies the second rule.
    * @details It eliminates recursively redundant states.
    * @note A redundant state is a state which has equivalent (same) successors
    * @warning This function invalidates the 'representative' map. After calling this function, we should recalculate it.
    * @param N the state who has potential redundant children
    * @param mapper a mapping between deleted redundant states and their replacement
    * */
        void join(node_t&N,std::unordered_map<node_t*,node_t*> &mapper)
        {
            //If N has a left child
            if(N.left)
            {
                //If this child has been previously removed, just update the child with its replacement
                if(mapper.count(N.left))
                    N.left=mapper[N.left];
                else
                {
                    //Else apply a recursive call
                    join(*N.left,mapper);
                    //If the child is redundant, remove it
                    if (N.left->left && N.left->left == N.left->right) {
                        auto tmp = N.left;
                        mapper[N.left] = N.left->left;
                        N.left = N.left->left;
                        delete tmp;
                    }
                }
            }
            if(N.right)
            {
                if(mapper.count(N.right))
                    N.right=mapper[N.right];
                else
                {
                    join(*N.right,mapper);
                    if (N.right->left && N.right->left == N.right->right) {
                        auto tmp = N.right;
                        mapper[N.right] = N.right->left;
                        N.right = N.right->left;
                        delete tmp;
                    }
                }
            }
        }

        void recalculate_representative(node_t*N)
        {
            if(!N)
                return;
            representative[*N]=N;
            recalculate_representative(N->left);
            recalculate_representative(N->right);
        }

        node_t* create_left(std::bitset<n> &S,int depth)
        {
            S[permutation[depth]]=false;
            return build(S,depth+1);
        }

        node_t* create_right(std::bitset<n> &S,int depth)
        {
            S[permutation[depth]]=true;
            return build(S,depth+1);
        }

        node_t* build(std::bitset<n> &S,int depth=0)
        {
            if(depth==n)
                return new node_t(evaluate(S));
            return new node_t(permutation[depth],create_left(S,depth),create_right(S,depth));
        }
        std::uint64_t count_truths() const {
            return FormalSpecification::count_truths(root,permutation.size());
        }

    };
}


#endif //VFSF_DECISION_TREE_H

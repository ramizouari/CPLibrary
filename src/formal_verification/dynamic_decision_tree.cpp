//
// Created by ramizouari on 18/10/22.
//

#include "formal_verification/dynamic_decision_tree.h"


namespace FormalSpecification
{
    dynamic_decision_tree::dynamic_decision_tree(int _n): abstract_decision_tree(_n),n(_n),permutation(_n)
    {
        for(int i=0;i<n;i++)
            permutation[i]=i;
    }

    dynamic_decision_tree::dynamic_decision_tree(std::vector<int> _permutation): abstract_decision_tree(_permutation.size()),permutation(std::move(_permutation)),n(permutation.size())
    {
    }

    dynamic_decision_tree::~dynamic_decision_tree()
    {
        for(auto &[_,node]:representative)
            delete node;
    }

    bool dynamic_decision_tree::evaluate_tree(const std::vector<bool> &S) const
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


    bool dynamic_decision_tree::operator()(const std::vector<bool> &S) const
    {
        return evaluate(S);
    }

    void dynamic_decision_tree::build()
    {
        std::vector<bool> S(n);
        root=build(S);
        compress_join();
    }

    natural dynamic_decision_tree::states_count() const
    {
        return representative.size();
    }

    natural dynamic_decision_tree::expanded_states_count() const
    {
        return (1<<(n+1))-1;
    }


    void dynamic_decision_tree::reseed(natural seed)
    {
        node_t::g.seed(seed);
    }


    void dynamic_decision_tree::compress_join()
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
     * */
    void dynamic_decision_tree::compress(node_t &N)
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
    void dynamic_decision_tree::join()
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
    void dynamic_decision_tree::join(node_t&N,std::unordered_map<node_t*,node_t*> &mapper)
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

    void dynamic_decision_tree::recalculate_representative(node_t*N)
    {
        if(!N)
            return;
        representative[*N]=N;
        recalculate_representative(N->left);
        recalculate_representative(N->right);
    }

    node_t* dynamic_decision_tree::create_left(std::vector<bool> &S,int depth)
    {
        S[permutation[depth]]=false;
        return build(S,depth+1);
    }

    node_t* dynamic_decision_tree::create_right(std::vector<bool> &S,int depth)
    {
        S[permutation[depth]]=true;
        return build(S,depth+1);
    }

    node_t* dynamic_decision_tree::build(std::vector<bool> &S,int depth)
    {
        if(depth==n)
            return new node_t(evaluate(S));
        return new node_t(permutation[depth],create_left(S,depth),create_right(S,depth));
    }

    dynamic_decision_tree::dynamic_decision_tree(const dynamic_decision_tree&O): abstract_decision_tree(O.n) {
        root=node_t::clone(O.root);
        permutation=O.permutation;
        recalculate_representative(root);
    }

    dynamic_decision_tree &dynamic_decision_tree::invert()
    {
        auto it1=representative.find(TrueNode),it2=representative.find(FalseNode);
        if(it1!=representative.end())
            *it1->second=FalseNode;
        if(it2!=representative.end())
            *it2->second=TrueNode;
        return *this;
    }


    bool default_dynamic_decision_tree::evaluate(const std::vector<bool> &S) const {
        return evaluate_tree(S);
    }


    default_dynamic_decision_tree operator^(const dynamic_decision_tree &A, const dynamic_decision_tree &B) {
        return merge(A,B,std::bit_xor<bool>());
    }

    default_dynamic_decision_tree operator|(const dynamic_decision_tree &A, const dynamic_decision_tree &B) {
        return merge(A,B,std::logical_or<>());
    }

    default_dynamic_decision_tree operator&(const dynamic_decision_tree &A, const dynamic_decision_tree &B)
    {
        return merge(A,B,std::logical_and<>());
    }

    default_dynamic_decision_tree operator~(const dynamic_decision_tree &A)
    {
        default_dynamic_decision_tree B(A.permutation);
        B.root=FormalSpecification::node_t::clone(A.root);
        //B.compress_join();
        B.recalculate_representative(B.root);
        B.invert();
        return B;
    }


    table_based_decision_tree::table_based_decision_tree(int n, const std::vector<bool> &_T): dynamic_decision_tree(n),T(_T)
    {
        if(T.size()!= 1LL<<this->n)
            throw std::runtime_error("Table Size is not a power of n");
    }

    table_based_decision_tree::table_based_decision_tree(int n, std::vector<bool> &&_T): dynamic_decision_tree(n),T(std::move(_T))
    {
        if(T.size()!= 1LL<<this->n)
            throw std::runtime_error("Table Size is not a power of n");
    }

    table_based_decision_tree::table_based_decision_tree(const std::vector<int> &permutation,
                                                         const std::vector<bool> &_T): dynamic_decision_tree(permutation),T(_T)
    {
        if(T.size()!= 1LL<<this->n)
            throw std::runtime_error("Table Size is not a power of n");
    }

    table_based_decision_tree::table_based_decision_tree(std::vector<int> &&permutation, const std::vector<bool> &_T): dynamic_decision_tree(std::move(permutation)),T(_T)
    {
        if(T.size()!= 1LL<<this->n)
            throw std::runtime_error("Table Size is not a power of n");
    }

    bool table_based_decision_tree::evaluate(const std::vector<bool> &S) const {
        std::uint64_t R=0;
        for(auto s:S)
            R=(R<<1)|static_cast<std::uint64_t>(s);
        return T[R];
    }

    table_based_decision_tree::table_based_decision_tree(const std::vector<int> &permutation, std::vector<bool> &&_T): dynamic_decision_tree(permutation),T(std::move(_T))
    {
        if(T.size()!= 1LL<<this->n)
            throw std::runtime_error("Table Size is not a power of n");
    }

    table_based_decision_tree::table_based_decision_tree(std::vector<int> &&permutation, std::vector<bool> &&_T): dynamic_decision_tree(std::move(permutation)),T(std::move(_T))
    {
        if(T.size()!= 1LL<<this->n)
            throw std::runtime_error("Table Size is not a power of n");
    }

    lambda_based_decision_tree::lambda_based_decision_tree(int n, const std::function<bool(std::vector<bool>)> &_T): dynamic_decision_tree(n),T(_T)
    {

    }

    lambda_based_decision_tree::lambda_based_decision_tree(const std::vector<int> &permutation,
                                                           const std::function<bool(std::vector<bool>)> &_T): dynamic_decision_tree(permutation),T(_T) {

    }

    lambda_based_decision_tree::lambda_based_decision_tree(const std::vector<int> &permutation,
                                                           std::function<bool(std::vector<bool>)> &&_T): dynamic_decision_tree(permutation),T(std::move(_T)) {

    }

    lambda_based_decision_tree::lambda_based_decision_tree(std::vector<int> &&permutation,
                                                           const std::function<bool(std::vector<bool>)> &_T):dynamic_decision_tree(std::move(permutation)),T(_T) {

    }

    lambda_based_decision_tree::lambda_based_decision_tree(std::vector<int> &&permutation,
                                                           std::function<bool(std::vector<bool>)> &&_T):dynamic_decision_tree(std::move(permutation)),T(std::move(_T)) {

    }

    bool lambda_based_decision_tree::evaluate(const std::vector<bool> &S) const {
        return T(S);
    }

    lambda_based_decision_tree::lambda_based_decision_tree(int n, std::function<bool(std::vector<bool>)>&&_T): dynamic_decision_tree(n),T(std::move(_T))
    {
    }
}
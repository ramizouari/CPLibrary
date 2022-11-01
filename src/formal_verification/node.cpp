//
// Created by ramizouari on 18/10/22.
//
#include "formal_verification/node.h"
#include <unordered_map>
namespace FormalSpecification
{
    node_t::node_t(natural symbol,node_t*left,node_t*right):symbol(symbol),left(left),right(right)
    {
        hash= (a*(left?left->hash:0)+b*symbol+c*(right?right->hash:0))%M;
    }
    node_t::node_t(Truth T):node_t(static_cast<natural>(T)){}

    node_t::node_t(bool T):node_t(T?Truth::True:Truth::False){}
    bool node_t::operator==(const node_t&O) const
    {
        if(hash!=O.hash || symbol!=O.symbol)
            return false;
        else if(left && !O.left || !left && O.left || left && O.left && *left!=*O.left)
            return false;
        else if(right && !O.right || !right && O.right || right && O.right && *right!=*O.right)
            return false;
        return true;
    }

    bool node_t::operator!=(const node_t&O) const
    {
        return !(*this == O);
    }

    node_t* clone(node_t*X,std::unordered_map<node_t*,node_t*>& mapper)
    {
        if(mapper.contains(X))
            return mapper[X];
        else return mapper[X]=new node_t(X->symbol,clone(X->left,mapper),clone(X->right,mapper));
    }

    node_t *node_t::clone(node_t*X) {
        //return new node_t(X->symbol,clone(X->left),clone(X->right));
        std::unordered_map<node_t*,node_t*> mapper;
        mapper[nullptr]= nullptr;
        return FormalSpecification::clone(X,mapper);
    }

    bool node_t::is_leaf() {
        return !left && !right;
    }


    struct TruthCache
    {
        std::uint64_t subtree_counter,subtree_variables;
    };

    std::uint64_t count_truths_rec(node_t*T ,std::uint64_t variables,std::unordered_map<node_t*,TruthCache> &cache)
    {
        if(cache.count(T))
        {
            auto [C,V]=cache[T];
            return C<<(variables-V);
        }
        else if(T->symbol==Truth::True)
        {
            cache[T]={1,0};
            return 1LL<<variables;
        }
        else if(T->symbol==Truth::False)
        {
            cache[T]={0,0};
            return 0;
        }
        auto &[C,V]=cache[T];
        C=count_truths_rec(T->left,variables-1,cache) + count_truths_rec(T->right,variables-1,cache);
        V=std::max(cache[T->left].subtree_variables,cache[T->right].subtree_variables)+1;
        return C;
    }

    std::uint64_t count_truths(node_t* T,int variables_count)
    {
        std::unordered_map<node_t*,TruthCache> C;
        return count_truths_rec(T,variables_count,C);
    }
}
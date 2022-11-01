//
// Created by ramizouari on 18/10/22.
//

#ifndef VFSF_NODE_H
#define VFSF_NODE_H
#include <cstdint>
#include <random>
namespace FormalSpecification
{
    using natural=std::uint64_t;

    enum Truth : natural
    {
        False=0xFFFFFFFFFFFFFFFE,True
    };

    struct node_t
    {
        inline static std::mt19937_64 g{};
        inline static constexpr natural M=1e9+7;
        inline static std::uniform_int_distribution<natural> d{1,M-1};
        inline static const natural a=d(g),b=d(g),c=d(g);
        static node_t* clone(node_t *S);
        natural hash;
        natural symbol;
        node_t *left,*right;
        explicit node_t(natural symbol,node_t*left=nullptr,node_t*right=nullptr);
        explicit node_t(Truth T);
        explicit node_t(bool T);
        bool operator==(const node_t&O) const;

        bool operator!=(const node_t&O) const;
        bool is_leaf();
    };
    const node_t TrueNode{Truth::True},FalseNode{Truth::False};

    std::uint64_t count_truths(node_t* T,int variables_count);

    template<typename Container, typename BooleanOperator>
    node_t* merge(node_t *u,node_t *v,Container && O , BooleanOperator && B)
    {
        if(u->is_leaf())
        {
            if(v->is_leaf())
                return new node_t(B(u->symbol==Truth::True,v->symbol==Truth::True));

            return new node_t(v->symbol,
                              merge(u,v->left,std::forward<Container>(O),std::forward<BooleanOperator>(B)),
                              merge(u,v->right,std::forward<Container>(O),std::forward<BooleanOperator>(B)));
        }
        else if(v->is_leaf() || O[u->symbol] > O[v->symbol])
            return merge(v,u,std::forward<Container>(O),std::forward<BooleanOperator>(B));
        else if(O[u->symbol] == O[v->symbol])
            return new node_t(u->symbol,
                              merge(u->left,v->left,std::forward<Container>(O),std::forward<BooleanOperator>(B)),
                              merge(u->right,v->right,std::forward<Container>(O),std::forward<BooleanOperator>(B)));
        else return new node_t(u->symbol,
                               merge(u->left,v,std::forward<Container>(O),std::forward<BooleanOperator>(B)),
                               merge(u->right,v,std::forward<Container>(O),std::forward<BooleanOperator>(B)));

    }
}



template<>
struct std::hash<FormalSpecification::node_t>
{
    using natural=FormalSpecification::natural;
    using node_t=FormalSpecification::node_t;
    inline static constexpr std::hash<int> hasher{};

    natural operator()(const node_t&N) const noexcept
    {
        return N.hash;
    }
    natural operator()(node_t*N) const
    {
        if(!N)
            return 0;
        return N->hash;
    }
};

#endif //VFSF_NODE_H

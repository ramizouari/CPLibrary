//
// Created by ramizouari on 07/11/23.
//

#ifndef CPLIBRARY_ISOMORPHISM_H
#define CPLIBRARY_ISOMORPHISM_H
#include "tree.h"

namespace graph
{
    using char_deque=std::deque<char>;

    std::shared_ptr<char_deque> string_encode(const graph::Tree & T, int u)
    {
        std::vector<std::shared_ptr<char_deque>> X;
        const auto &C=T.children(u);
        if(C.empty())
            return std::make_shared<char_deque>(char_deque{'(', ')'});
        X.reserve(C.size());
        for(auto v:C) X.push_back(string_encode(T,v));
        std::sort(X.begin(),X.end(),[](const auto &x,const auto &y)
        {
            return *x < *y;
        });
        auto it=std::max_element(X.begin(),X.end(),[](const auto &x,const auto &y) {
            return x->size() < y->size();
        });
        auto r=std::distance(X.begin(),it);
        std::shared_ptr<char_deque> Z=X[r];
        for(int i=r-1;i>=0;i--)
            std::copy(X[i]->rbegin(),X[i]->rend(),std::front_inserter(*Z));
        for(int i=r+1;i<X.size();i++)
            std::copy(X[i]->begin(),X[i]->end(),std::back_inserter(*Z));
        Z->push_back(')');
        Z->push_front('(');
        return Z;
    }

    std::string string_encode(const graph::Tree &T)
    {
        auto E=string_encode(T,T.root);
        return std::string(E->begin(),E->end());
    }

    std::optional<int> second_centroid(graph::Tree & T)
    {
        auto u=T.root;
        auto X=T.children(u);
        for(auto v:X)
        {
            bool is_centroid=true;
            T.adjacentReRoot(v);
            if(T.subtree_size[u]>T.subtree_size[v]/2)
                is_centroid=false;
            T.adjacentReRoot(u);
            if(is_centroid)
                return v;
        }
        return std::nullopt;
    }

    std::vector<std::string> full_string_encoding(graph::Tree &T)
    {
        T.centroid();
        T.buildStatistics(graph::TreeStats::SIZE);
        std::vector<std::string> A;
        A.push_back(string_encode(T));
        auto p= second_centroid(T);
        if(p.has_value())
        {
            T.reRoot(*p);
            T.buildStatistics(graph::TreeStats::SIZE);
            A.push_back(string_encode(T));
            if(A[0] > A[1])
                std::swap(A[0],A[1]);
        }
        return A;
    }

    struct TreeHolder
    {
        std::shared_ptr<Tree> tree;
        TreeHolder(Tree&& t): tree(std::make_shared<Tree>(std::move(t))){}
        TreeHolder(int n): tree(std::make_shared<Tree>(n)){}

        Tree* operator->() const
        {
            return tree.get();
        }
        Tree& operator*() const
        {
            return *tree;
        }
    };


    struct IsoTreeEq
    {
        bool operator()(Tree &a, Tree &b) const
        {
            if(a.size() != b.size())
                return false;
            a.centroid();
            b.centroid();
            auto M1= string_encode(a);
            auto M2= string_encode(b);
            if(M1==M2)
                return true;
            auto c= second_centroid(a);
            if(c.has_value())
            {
                a.reRoot(*c);
                a.buildStatistics(graph::TreeStats::SIZE);
                auto M3= string_encode(a);
                return M3==M2;
            }
            return false;
        }

        bool operator()(const TreeHolder &a, const TreeHolder &b) const
        {
            return (*this)(*a,*b);
        }
    };

    struct IsoTreeCmp
    {
        bool operator()(Tree &a, Tree &b) const
        {
            if(a.size()!=b.size())
                return a.size() < b.size();
            auto X=full_string_encoding(a);
            auto Y=full_string_encoding(b);
            for(int i:{0,1}) for(int j:{0,1}) if(i<X.size() && j<Y.size() && X[i]==Y[j])
            {
                if(i!=0 || j!=0)
                    throw std::runtime_error("?????????????????????????????????????");
                else return false;
            }
            return X < Y;
        }

        bool operator()(const TreeHolder& a, const TreeHolder& b) const
        {
            return (*this)(*a,*b);
        }
    };

    struct IsoTreeHash
    {
        std::hash<std::string> H;
    public:
        size_t operator()(Tree &a) const
        {
            auto A= full_string_encoding(a);
            return std::accumulate(A.begin(),A.end(),0ULL,[&H=this->H](auto x,auto y){
                return x^H(y);
            });
        }

        size_t operator()(const TreeHolder& a) const
        {
            return (*this)(*a);
        }
    };
}

#endif //CPLIBRARY_ISOMORPHISM_H

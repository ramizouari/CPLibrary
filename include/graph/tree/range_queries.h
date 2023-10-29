//
// Created by ramizouari on 25/10/23.
//

#ifndef CPLIBRARY_RANGE_QUERIES_H
#define CPLIBRARY_RANGE_QUERIES_H
#include "tree.h"

namespace graph
{
    template<typename O, typename RMQ>
    struct HeavyLightTree : public WeightedTree<typename O::type>
    {
        using type=typename O::type;
        using Weight=type;
    private:
    public:
        using WeightedTree<Weight>::HLD;
        using WeightedTree<Weight>::WeightedTree;
        using AdjacentType=WeightedTree<Weight>::AdjacentType;
        using WeightedTree<Weight>::root;
        using WeightedTree<Weight>::parent;
        inline static O F{};

        void buildStatistics(TreeStats stats=TreeStats::RANGE_QUERIES)
        {
            WeightedTree<Weight>::buildStatistics(stats);
            if(stats & TreeStats::RANGE_QUERIES) buildRangeStatistics();
        }

        Weight query(int u,int v)
        {
            auto lca=WeightedTree<Weight>::leastCommonAncestor(u,v);
            return F(query_with_lca(u,lca,false),query_with_lca(v,lca,true));
        }

        void update(int u,int v, Weight w)
        {
            if(parent[u] && parent[u]->first==v)
                return update(u,w);
            if(parent[v] && parent[v]->first==u)
                return update(v,w);
            throw std::invalid_argument("u and v must be adjacent");
        }

        void update(int u,Weight w)
        {
            if(!parent[u])
                throw std::invalid_argument("u must not be the root");
            parent[u]->second=w;
            auto v=parent[u]->first;
            if(WeightedTree<Weight>::HLD.is_heavy[u])
            {
                auto &R_left=S_left[HLD.HLD_mapper[v].hld_id];
                auto &R_right=S_right[HLD.HLD_mapper[v].hld_id];
                R_left.update(HLD.component_size[HLD.HLD_mapper[v].hld_id] - HLD.HLD_mapper[v].index - 1, w);
                R_right.update(HLD.HLD_mapper[v].index, w);
            }
        }

    protected:


        //TODO: Support for non-commutative operations
        Weight query_with_lca(int u,int lca, bool invert)
        {
            int o=u;
            if(u==lca)
                return O::neutral;
            Weight R=O::neutral;
            while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id)
            {
                auto [x,y]=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id];
                auto b=HLD.HLD_mapper[u].index;
                auto a= HLD.HLD_mapper[x].index;
                auto &S=invert?S_right:S_left;
                if(invert)
                    R=F(S[HLD.HLD_mapper[u].hld_id].query(a, b),R);
                else
                    R=F(R,S[HLD.HLD_mapper[u].hld_id].query(a, b));
                u=x;
                while(u!=lca && !WeightedTree<Weight>::HLD.is_heavy[u])
                {
                    if(invert) R=F(parent[u]->second,R);
                    else R=F(R,parent[u]->second);
                    u=parent[u]->first;
                }
            }
            auto b=HLD.HLD_mapper[u].index;
            auto a= HLD.HLD_mapper[lca].index;
            auto &S=invert?S_right:S_left;
            if(invert)
                R=F(S[HLD.HLD_mapper[u].hld_id].query(a, b),R);
            else
                R=F(R,S[HLD.HLD_mapper[u].hld_id].query(a, b));
            return R;
        }

        std::vector<RMQ> S_left, S_right;

        void buildRangeStatistics()
        {
            for(auto &C:HLD.components)
            {
                HLD.component_size.push_back(C.size());
                if(C.empty()) C.emplace_back(O::neutral);
                S_right.emplace_back(C);
                std::reverse(C.begin(),C.end());
                S_left.emplace_back(std::move(C));
            }
        }
    };

    template<typename O, typename RMQ>
    struct HeavyLightNodeTree : protected HeavyLightTree<O,RMQ>
    {
        using type=typename O::type;
        using Weight=type;
    protected:
        void buildRepresentation(int root)
        {
            reRoot(root);
            for(int u=0;u<n;u++) for(auto &[v,w]: this->reverseList[u])
                {
                    w=node_weight[v];
                    parent[v]=std::make_pair(u,w);
                    this->adjacencyList[v].front().second=w;
                }
        }
        using WeightedTree<Weight>::children;
        using HeavyLightTree<O,RMQ>::query_with_lca;
        using HeavyLightTree<O,RMQ>::F;
        using WeightedTree<Weight>::reRoot;
    public:
        std::vector<Weight> node_weight;
        HeavyLightNodeTree(int n,int root):HeavyLightTree<O,RMQ>(n,root),node_weight(n,O::neutral){}
        HeavyLightNodeTree(int n):HeavyLightNodeTree(n,0){}
        using AdjacentType=WeightedTree<Weight>::AdjacentType;
        using WeightedTree<Weight>::root;
        using WeightedTree<Weight>::parent;
        using WeightedTree<Weight>::n;

        void setWeight(int u,Weight w)
        {
            node_weight[u]=w;
        }

        void connect(int u,int v)
        {
            HeavyLightTree<O,RMQ>::connect(u,v,O::neutral);
        }

        void setParent(int u,int v)
        {
            HeavyLightTree<O,RMQ>::setParent(u,v,O::neutral);
        }

        void build(int root=0,TreeStats stats= TreeStats::RANGE_QUERIES)
        {
            buildRepresentation(root);
            HeavyLightTree<O,RMQ>::buildStatistics(stats);
        }

        Weight query(int u,int v)
        {
            auto lca=WeightedTree<Weight>::leastCommonAncestor(u,v);
            return F(query_with_lca(u,lca,false),node_weight[lca],query_with_lca(v,lca,true));
        }

        void update(int u, Weight w)
        {
            node_weight[u]=w;
            if(parent[u])
                HeavyLightTree<O, RMQ>::update(u, w);
        }

        using HeavyLightTree<O,RMQ>::distance;
    };


    template<typename O, typename RMQ>
    struct CommutativeHeavyLightTree : public WeightedTree<typename O::type>
    {
        using Weight=typename O::type;
        using type=Weight;
    private:
    public:
        using WeightedTree<Weight>::HLD;
        using WeightedTree<Weight>::WeightedTree;
        using AdjacentType=WeightedTree<Weight>::AdjacentType;
        using WeightedTree<Weight>::root;
        using WeightedTree<Weight>::parent;
        inline static O F{};

        void buildStatistics(TreeStats stats=TreeStats::RANGE_QUERIES)
        {
            WeightedTree<Weight>::buildStatistics(stats);
            if(stats & TreeStats::RANGE_QUERIES) buildRangeStatistics();
        }

        Weight query(int u,int v)
        {
            auto lca=WeightedTree<Weight>::leastCommonAncestor(u,v);
            return F(query_with_lca(u,lca),query_with_lca(v,lca));
        }

        void update(int u,int v, Weight w)
        {
            if(parent[u] && parent[u]->first==v)
                return update(u,w);
            if(parent[v] && parent[v]->first==u)
                return update(v,w);
            throw std::invalid_argument("u and v must be adjacent");
        }

        void update(int u,Weight w)
        {
            if(!parent[u])
                throw std::invalid_argument("u must not be the root");
            parent[u]->second=w;
            auto v=parent[u]->first;
            if(WeightedTree<Weight>::HLD.is_heavy[u])
            {
                auto &R=S[HLD.HLD_mapper[v].hld_id];
                R.update(HLD.HLD_mapper[v].index, w);
            }
        }

    protected:
        Weight query_with_lca(int u,int lca)
        {
            if(u==lca)
                return O::neutral;
            Weight R=O::neutral;
            while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id)
            {
                auto [x,y]=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id];
                auto b=HLD.HLD_mapper[u].index;
                auto a= HLD.HLD_mapper[x].index;
                R=F(S[HLD.HLD_mapper[u].hld_id].query(a, b), R);
                u=x;
                while(u!=lca && !WeightedTree<Weight>::HLD.is_heavy[u])
                {
                    R=F(R,parent[u]->second);
                    u=parent[u]->first;
                }
            }
            auto b=HLD.HLD_mapper[u].index;
            auto a= HLD.HLD_mapper[lca].index;
            R=F(R,S[HLD.HLD_mapper[u].hld_id].query(a, b));
            return R;
        }

        std::vector<RMQ> S;

        void buildRangeStatistics()
        {
            for(auto &C:HLD.components)
            {
                HLD.component_size.push_back(C.size());
                if(C.empty()) C.emplace_back(O::neutral);
                S.emplace_back(C);
            }
        }
    };



    template<typename O, typename RMQ>
    struct CommutativeHeavyLightNodeTree : protected CommutativeHeavyLightTree<O,RMQ>
    {
        using type=typename O::type;
        using Weight=type;
    protected:
        void buildRepresentation(int root)
        {
            reRoot(root);
            for(int u=0;u<n;u++) for(auto &[v,w]: this->reverseList[u])
                {
                    w=node_weight[v];
                    parent[v]=std::make_pair(u,w);
                    this->adjacencyList[v].front().second=w;
                }
        }
        using WeightedTree<Weight>::children;
        using CommutativeHeavyLightTree<O,RMQ>::query_with_lca;
        using CommutativeHeavyLightTree<O,RMQ>::F;
        using WeightedTree<Weight>::reRoot;
    public:
        std::vector<Weight> node_weight;
        CommutativeHeavyLightNodeTree(int n,int root): CommutativeHeavyLightTree<O,RMQ>(n,root),node_weight(n,O::neutral){}
        CommutativeHeavyLightNodeTree(int n):CommutativeHeavyLightNodeTree(n,0){}
        using AdjacentType=WeightedTree<Weight>::AdjacentType;
        using WeightedTree<Weight>::root;
        using WeightedTree<Weight>::parent;
        using WeightedTree<Weight>::n;

        void setWeight(int u,Weight w)
        {
            node_weight[u]=w;
        }

        void connect(int u,int v)
        {
            CommutativeHeavyLightTree<O,RMQ>::connect(u,v,O::neutral);
        }

        void setParent(int u,int v)
        {
            CommutativeHeavyLightTree<O,RMQ>::setParent(u,v,O::neutral);
        }

        void build(int root=0,TreeStats stats= TreeStats::RANGE_QUERIES)
        {
            buildRepresentation(root);
            CommutativeHeavyLightTree<O,RMQ>::buildStatistics(stats);
        }

        Weight query(int u,int v)
        {
            auto lca=WeightedTree<Weight>::leastCommonAncestor(u,v);
            return F(query_with_lca(u,lca),node_weight[lca],query_with_lca(v,lca));
        }

        void update(int u, Weight w)
        {
            node_weight[u]=w;
            if(parent[u])
                CommutativeHeavyLightTree<O, RMQ>::update(u, w);
        }

        using CommutativeHeavyLightTree<O,RMQ>::distance;
    };
}

#endif //CPLIBRARY_RANGE_QUERIES_H

//
// Created by ramizouari on 24/10/23.
//

#ifndef CPLIBRARY_TREE_H
#define CPLIBRARY_TREE_H

#include <optional>
#include <stdexcept>
#include "graph/graph.h"
#include "data_structures/dynamic/sparse_array.h"
#include "data_structures/fixed/sparse_array.h"
#include "algebra/binary_operation.h"
#include <memory>
#include <queue>
#include <random>
#include <chrono>

namespace cp::graph
{

    namespace ds= data_structures;

    struct HLDIndex
    {
        int hld_id;
        int index;
        HLDIndex(int _hld_id, int _index): hld_id(_hld_id), index(_index){}
    };

    template<typename Weight>
    struct HeavyLightDecomposition
    {
        std::vector<bool> is_heavy;
        std::vector<std::pair<int,int>> heavy_path_endpoints;
        std::vector<int> component_size;
        std::vector<HLDIndex> HLD_mapper;
        std::vector<std::vector<Weight>> components;
    };

    template<>
    struct HeavyLightDecomposition<void>
    {
        std::vector<bool> is_heavy;
        std::vector<std::pair<int,int>> heavy_path_endpoints;
        std::vector<int> component_size;
        std::vector<HLDIndex> HLD_mapper;
    };

    enum class TreeStats
    {
        NONE=0,
        SIZE=0b00001,
        HEAVY_EDGES=0b00011,
        LCA=0b00100,
        HLD=0b01111,
        RANGE_QUERIES=0b11111,
    };
    bool operator &(TreeStats a,TreeStats b)
    {
        return static_cast<int>(a) & static_cast<int>(b);
    }

    struct Tree : Graph
    {

        bool reversed=false;
        std::vector<int> subtree_size;
        std::vector<std::optional<int>> parent;
        int root;
        HeavyLightDecomposition<void> HLD;
        Tree(int n,int _root):Graph(n),root(_root),subtree_size(n),parent(n)
        {
            HLD.is_heavy.resize(n,false);
        }

        Tree(const Graph &G): Graph(G),root(G.n),subtree_size(G.n),parent(G.n)
        {
            HLD.is_heavy.resize(G.n,false);
        }

        Tree(Graph &&G): Graph(std::move(G)),root(n),subtree_size(n),parent(n)
        {
            HLD.is_heavy.resize(n,false);
        }

        explicit Tree(int n):Tree(n,0){}

        void setParent(int u,int v)
        {
            if(reversed) std::swap(u,v);
            this->connect(u,v);
            parent[u].emplace(v);
        }

        std::vector<int> &children(int u)
        {
            return reversed?this->adjacencyList[u] : this->reverseList[u] ;
        }

        [[nodiscard]] const std::vector<int> &children(int u) const
        {
            return reversed?this->adjacencyList[u] : this->reverseList[u] ;
        }

        void buildStatistics(TreeStats stats=TreeStats::HLD)
        {
            updateRoot();
            if(stats & TreeStats::SIZE) updateSize(root);
            //TODO: Optimize heavy edges
            if(stats & TreeStats::HEAVY_EDGES) updateHeavyEdges(root);
            if(stats & TreeStats::LCA) buildLCA();
            auto t1=std::chrono::high_resolution_clock::now();
            //TODO: Optimize HLD
            if(stats & TreeStats::HLD) buildHeavyLightDecomposition();
        }

        void adjacentReRoot(int new_root)
        {
            if(parent[new_root] != root)
                throw std::invalid_argument("new root must be adjacent to old root");
            auto u=*parent[new_root];
            parent[new_root]=std::nullopt;
            parent[root].emplace(new_root);
            auto delta=subtree_size[new_root];
            subtree_size[new_root]=subtree_size[root];
            subtree_size[root]-=delta;
            root=new_root;
        }

        void reRoot(int new_root)
        {
            std::queue<int> Q;
            std::vector<bool> visited(n);
            Q.push(new_root);
            visited[new_root]=true;
            std::vector<std::vector<int>> newAdjacencyList(n),newReverseList(n);
            while(!Q.empty())
            {
                auto u=Q.front();
                Q.pop();
                for(auto v:this->adjacencyList[u]) if(!visited[v])
                    {
                        visited[v]=true;
                        Q.push(v);
                        newReverseList[u].emplace_back(v);
                        newAdjacencyList[v].emplace_back(u);
                    }
                for(auto v:this->reverseList[u]) if(!visited[v])
                    {
                        visited[v]=true;
                        Q.push(v);
                        newReverseList[u].emplace_back(v);
                        newAdjacencyList[v].emplace_back(u);
                    }
            }
            this->adjacencyList=std::move(newAdjacencyList);
            this->reverseList=std::move(newReverseList);
            updateRoot();
        }

        int leastCommonAncestor(int u,int v)
        {
            if(lca_data)
            {
                auto [a,b]=euler_tour_endpoints[u];
                auto [c,d]=euler_tour_endpoints[v];
                if(a>c)
                {
                    std::swap(a,c);
                    std::swap(b,d);
                }
                if(b<d)
                    return lca_data->query(a,c).second;
                else
                    return lca_data->query(a,d).second;
            }
            else
            {
                while(u!=v)
                {
                    if(subtree_size[u]>subtree_size[v])
                        u=*parent[u];
                    else
                        v=*parent[v];
                }
                return u;
            }
        }
        void updateSize(int u)
        {
            subtree_size[u]=1;
            for(auto v:children(u))
            {
                updateSize(v);
                subtree_size[u]+=subtree_size[v];
            }
        }

        void updateRoot()
        {
            for(int i=0;i<n;i++)
                parent[i]=this->adjacencyList[i].empty()?std::nullopt:std::make_optional(this->adjacencyList[i][0]);
            for(int i=0;i< n;i++)
                if(!parent[i])
                    root=i;
        }

        void updateHeavyEdges(int u)
        {
            for(int i=0;i<children(u).size();i++)
            {
                auto v=children(u)[i];
                if(subtree_size[v]>=(subtree_size[u]+1)/2)
                    HLD.is_heavy[v]=true;
                updateHeavyEdges(v);
            }
        }
        void buildLCA()
        {
            std::vector<HeightData> A;
            euler_tour_endpoints.resize(n);
            eulerTour(root,0,A);
            min_t<HeightData>::neutral.first=std::numeric_limits<int>::max();
            lca_data=std::make_unique<ds::fixed::sparse_array<min_t<HeightData>>>(A);
        }

        void buildHeavyLightDecomposition()
        {
            std::vector<int> stack;
            stack.push_back(root);
            HLD.HLD_mapper.resize(n, HLDIndex(-1, -1));
            HLD.HLD_mapper[root]={0, 0};
            HLD.heavy_path_endpoints.emplace_back(root,root);
            int components=1;
            while(!stack.empty())
            {
                auto u=stack.back();
                stack.pop_back();
                for(auto v:children(u))
                {
                    if(HLD.is_heavy[v])
                    {
                        auto &[_,y] = HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id];
                        stack.push_back(v);
                        HLD.HLD_mapper[v]={HLD.HLD_mapper[u].hld_id, HLD.HLD_mapper[u].index + 1};
                        y=v;
                    }
                    else
                    {
                        HLD.heavy_path_endpoints.emplace_back(v,v);
                        stack.push_back(v);
                        HLD.HLD_mapper[v]=HLDIndex(components, 0);
                        components++;
                    }
                }
            }
        }

        int distance(int a,int b)
        {
            auto lca=leastCommonAncestor(a,b);
            return distance_with_lca(a,lca)+distance_with_lca(b,lca);
        }

        int distance_with_lca(int u,int lca)
        {
            int d=0;
            while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id)
            {
                if(HLD.is_heavy[u])
                {
                    d+=HLD.HLD_mapper[u].index;
                    u=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id].first;
                }
                else
                {
                    d++;
                    u=*parent[u];
                }
            }
            return d+= HLD.HLD_mapper[u].index - HLD.HLD_mapper[lca].index;
        }

        int centroid()
        {
            reRoot(centroid(root,root));
            return root;
        }

        int centroid(int u)
        {
            return centroid(u,std::nullopt);
        }

        int centroid(int u,std::optional<int> p)
        {
            for(auto v:children(u)) if(v!= p && subtree_size[v]>=(subtree_size[u]+1)/2)
                {
                    adjacentReRoot(v);
                    return centroid(v,u);
                }
            return u;
        }
    protected:
        using HeightData=std::pair<int,int>;
        using EnpointsData = std::pair<int,int>;
        std::unique_ptr<ds::fixed::sparse_array<min_t<HeightData>>> lca_data;
        std::vector<EnpointsData> euler_tour_endpoints;
        void eulerTour(int u,int height,std::vector<HeightData> &A)
        {
            euler_tour_endpoints[u].first=A.size();
            for(auto v: children(u))
            {
                A.emplace_back(height,u);
                eulerTour(v,height+1,A);
            }
            A.emplace_back(height,u);
            euler_tour_endpoints[u].second=A.size();
        }
    };

    template<typename Weight>
    struct WeightedTree : public WeightedGraph<Weight>
    {

        bool reversed=false;
        std::vector<int> subtree_size;
        using AdjacentType=WeightedGraph<Weight>::AdjacentType;
        std::vector<std::optional<AdjacentType>> parent;
        int root;
        HeavyLightDecomposition<Weight> HLD;
        WeightedTree(int n,int _root):WeightedGraph<Weight>(n),root(_root),subtree_size(n),parent(n)
        {
            HLD.is_heavy.resize(n);
        }
        explicit WeightedTree(int n):WeightedTree(n,0){}

        void setParent(int u,int v,Weight w)
        {
            if(reversed) std::swap(u,v);
            this->connect(u,v,w);
            parent[u].emplace(v,w);
        }

        std::vector<AdjacentType> &children(int u)
        {
            return reversed?this->adjacencyList[u] : this->reverseList[u] ;
        }

        void buildStatistics(TreeStats stats=TreeStats::HLD)
        {
            updateRoot();
            if(stats & TreeStats::SIZE) updateSize(root);
            if(stats & TreeStats::HEAVY_EDGES) updateHeavyEdges(root);
            if(stats & TreeStats::LCA) buildLCA();
            if(stats & TreeStats::HLD) buildHeavyLightDecomposition();
        }

        void adjacentReRoot(int new_root)
        {
            if(!parent[new_root] || parent[new_root]->first != root)
                throw std::invalid_argument("new root must be adjacent to old root");
            auto [u,w]=*parent[new_root];
            parent[new_root]=std::nullopt;
            parent[root].emplace(new_root,w);
            auto delta=subtree_size[new_root];
            subtree_size[new_root]=subtree_size[root];
            subtree_size[root]-=delta;
        }

        void reRoot(int new_root)
        {
            std::queue<int> Q;
            std::vector<bool> visited(WeightedGraph<Weight>::n);
            Q.emplace(new_root);
            visited[new_root]=true;
            std::vector<std::vector<AdjacentType>> newAdjacencyList(WeightedGraph<Weight>::n),newReverseList(WeightedGraph<Weight>::n);
            while(!Q.empty())
            {
                auto u=Q.front();
                Q.pop();
                for(auto [v,w]:this->adjacencyList[u]) if(!visited[v])
                    {
                        visited[v]=true;
                        Q.emplace(v);
                        newReverseList[u].emplace_back(v,w);
                        newAdjacencyList[v].emplace_back(u,w);
                    }
                for(auto [v,w]:this->reverseList[u]) if(!visited[v])
                    {
                        visited[v]=true;
                        Q.emplace(v);
                        newReverseList[u].emplace_back(v,w);
                        newAdjacencyList[v].emplace_back(u,w);
                    }
            }
            this->adjacencyList=std::move(newAdjacencyList);
            this->reverseList=std::move(newReverseList);
            updateRoot();
        }

        int leastCommonAncestor(int u,int v)
        {
            if(lca_data)
            {
                auto [a,b]=euler_tour_endpoints[u];
                auto [c,d]=euler_tour_endpoints[v];
                if(a>c)
                {
                    std::swap(a,c);
                    std::swap(b,d);
                }
                if(b<d)
                    return lca_data->query(a,c).second;
                else
                    return lca_data->query(a,d).second;
            }
            else
            {
                while(u!=v)
                {
                    if(subtree_size[u]>subtree_size[v])
                        u=parent[u]->first;
                    else
                        v=parent[v]->first;
                }
                return u;
            }
        }
        void updateSize(int u)
        {
            subtree_size[u]=1;
            for(auto [v,_]:children(u))
            {
                updateSize(v);
                subtree_size[u]+=subtree_size[v];
            }
        }

        void updateRoot()
        {
            for(int i=0;i<WeightedGraph<Weight>::n;i++)
                parent[i]=this->adjacencyList[i].empty()?std::nullopt:std::make_optional(this->adjacencyList[i][0]);
            for(int i=0;i<WeightedGraph<Weight>::n;i++)
                if(!parent[i])
                    root=i;
        }

        void updateHeavyEdges(int u)
        {
            for(int i=0;i<children(u).size();i++)
            {
                auto [v,_]=children(u)[i];
                if(subtree_size[v]>=(subtree_size[u]+1)/2)
                    HLD.is_heavy[v]=true;
                updateHeavyEdges(v);
            }
        }
        void buildLCA()
        {
            std::vector<HeightData> A;
            euler_tour_endpoints.resize(WeightedGraph<Weight>::n);
            eulerTour(root,0,A);
            min_t<HeightData>::neutral.first=std::numeric_limits<int>::max();
            lca_data=std::make_unique<ds::fixed::sparse_array<min_t<HeightData>>>(A);
        }

        void buildHeavyLightDecomposition()
        {
            std::vector<int> stack;
            stack.push_back(root);
            HLD.HLD_mapper.resize(WeightedTree<Weight>::n, HLDIndex(-1, -1));
            HLD.HLD_mapper[root]={0, 0};
            HLD.components.emplace_back();
            HLD.heavy_path_endpoints.emplace_back(root,root);
            while(!stack.empty())
            {
                auto u=stack.back();
                stack.pop_back();
                for(auto [v,w]:children(u))
                {
                    if(HLD.is_heavy[v])
                    {
                        auto &C=HLD.components[HLD.HLD_mapper[u].hld_id];
                        auto &[_,y] = HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id];
                        stack.push_back(v);
                        HLD.HLD_mapper[v]={HLD.HLD_mapper[u].hld_id, HLD.HLD_mapper[u].index + 1};
                        y=v;
                        C.push_back(w);
                    }
                    else
                    {
                        HLD.heavy_path_endpoints.emplace_back(v,v);
                        stack.push_back(v);
                        HLD.HLD_mapper[v]=HLDIndex(HLD.components.size(), 0);
                        HLD.components.emplace_back();
                    }
                }
            }
            for(const auto &C:HLD.components)
                HLD.component_size.push_back(C.size());
        }

        int distance(int a,int b)
        {
            auto lca=leastCommonAncestor(a,b);
            return distance_with_lca(a,lca)+distance_with_lca(b,lca);
        }

        int distance_with_lca(int u,int lca)
        {
            int d=0;
            while(HLD.HLD_mapper[u].hld_id != HLD.HLD_mapper[lca].hld_id)
            {
                if(WeightedTree<Weight>::HLD.is_heavy[u])
                {
                    d+=HLD.HLD_mapper[u].index;
                    u=HLD.heavy_path_endpoints[HLD.HLD_mapper[u].hld_id].first;
                }
                else
                {
                    d++;
                    u=parent[u]->first;
                }
            }
            return d+= HLD.HLD_mapper[u].index - HLD.HLD_mapper[lca].index;
        }


    protected:
        using HeightData=std::pair<int,int>;
        using EnpointsData = std::pair<int,int>;
        std::unique_ptr<ds::fixed::sparse_array<min_t<HeightData>>> lca_data;
        std::vector<EnpointsData> euler_tour_endpoints;
        void eulerTour(int u,int height,std::vector<HeightData> &A)
        {
            euler_tour_endpoints[u].first=A.size();
            for(auto [v,_]: children(u))
            {
                A.emplace_back(height,u);
                eulerTour(v,height+1,A);
            }
            A.emplace_back(height,u);
            euler_tour_endpoints[u].second=A.size();
        }
    };

}


#endif //CPLIBRARY_TREE_H

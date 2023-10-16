//
// Created by ramizouari on 21/11/22.
//

#ifndef CPLIBRARY_UNION_FIND_H
#define CPLIBRARY_UNION_FIND_H
#include <vector>
#include <map>
#include <unordered_map>

class UnionFind
{
    int n;
    std::vector<int> parent,rank;
public:
    UnionFind(int n):n(n),rank(n),parent(n)
    {
        for(int i=0;i<n;i++)
            parent[i]=i;
    }

    void connect(int a,int b)
    {
        auto u= representative(a),v= representative(b);
        if(u==v)
            return;
        if(rank[u]<rank[u])
            parent[u]=parent[v];
        else if(rank[v]<rank[u])
            parent[v]=parent[u];
        else
        {
            parent[u]=parent[v];
            rank[v]++;
        }
    }

    int representative(int a)
    {
        if(parent[a]==a)
            return a;
        else return parent[a]= representative(parent[a]);
    }

    bool equivalent(int a,int b)
    {
        return representative(a)== representative(b);
    }
};

template<typename OrderedSet>
class OrderedUnionFind
{
    struct ParentData
    {
        const OrderedSet* parent;
        int rank;
    };
    std::map<OrderedSet,ParentData> data{};
public:
    void connect(const OrderedSet& A,const OrderedSet&B)
    {
        if(equivalent(A,B))
            return;
        auto &[C,NodeData1]=get(representative(A));
        auto &[D,NodeData2]=get(representative(B));
        auto &[p1,r1]=NodeData1;
        auto &[p2,r2]=NodeData2;
        if(r1<r2)
            p1=p2;
        else if(r1>r2)
            p2=p1;
        else
        {
            p1=p2;
            r2++;
        }
    }
    std::pair<const OrderedSet,ParentData>& get(const OrderedSet &A)
    {
        if(!data.contains(A))
        {
            auto [it,_]=data.emplace(A,ParentData{nullptr,0});
            auto &[B,NodeData]=*it;
            NodeData.parent=&B;
            return *it;
        }
        return *data.find(A);
    }
    const OrderedSet& representative(const OrderedSet&A)
    {
        auto &[B,NodeData]=get(A);
        auto &C=NodeData.parent;
        if(&B==C)
            return B;
        else {
            NodeData.parent = &representative(*NodeData.parent);
            return *NodeData.parent;
        }
    }

    bool equivalent(const OrderedSet&A,const OrderedSet&B)
    {
        return &representative(get(A).first)== &representative(get(B).first);
    }

};

template<typename UnorderedSet>
class UnorderedUnionFind
{
    struct ParentData
    {
        const UnorderedSet* parent;
        int rank;
    };
    std::unordered_map<UnorderedSet,ParentData> data{};
public:
    void connect(const UnorderedSet& A,const UnorderedSet&B)
    {
        if(equivalent(A,B))
            return;
        auto &[C,NodeData1]=get(representative(A));
        auto &[D,NodeData2]=get(representative(B));
        auto &[p1,r1]=NodeData1;
        auto &[p2,r2]=NodeData2;
        if(r1<r2)
            p1=p2;
        else if(r1>r2)
            p2=p1;
        else
        {
            p1=p2;
            r2++;
        }
    }
    std::pair<const UnorderedSet,ParentData>& get(const UnorderedSet &A)
    {
        if(!data.contains(A))
        {
            auto [it,_]=data.emplace(A,ParentData{nullptr,0});
            auto &[B,NodeData]=*it;
            NodeData.parent=&B;
            return *it;
        }
        return *data.find(A);
    }
    const UnorderedSet& representative(const UnorderedSet&A)
    {
        auto &[B,NodeData]=get(A);
        auto &C=NodeData.parent;
        if(&B==C)
            return B;
        else {
            NodeData.parent = &representative(*NodeData.parent);
            return *NodeData.parent;
        }
    }

    bool equivalent(const UnorderedSet&A,const UnorderedSet&B)
    {
        return &representative(get(A).first)== &representative(get(B).first);
    }

};



#endif //CPLIBRARY_UNION_FIND_H

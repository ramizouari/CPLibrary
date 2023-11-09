#include "graph/tree/tree.h"
#include "graph/tree/isomorphism.h"
#include <iostream>

using TreeesHashMap=std::unordered_map<graph::TreeHolder,double,graph::IsoTreeHash,graph::IsoTreeEq>;

struct NaiveTreeCmp
{
    bool operator()(const graph::TreeHolder &a, const graph::TreeHolder &b) const
    {
        return a->adjacencyList < b->adjacencyList;
    }
};
using TreesMap=std::map<graph::TreeHolder,double,graph::IsoTreeCmp>;
using NaiveMap=std::map<graph::TreeHolder,double,NaiveTreeCmp>;


graph::TreeHolder cut(const graph::Tree & A,int node)
{
    auto n=A.size();
    graph::Tree B(n-1,0);
    for(int i=0;i<n;i++) if(i!=node) for(auto j:A.children(i)) if(j!=node)
        B.setParent(j-(j>node),i-(i>node));
    B.reRoot(0);
    B.buildStatistics(graph::TreeStats::SIZE);
    return B;
}

std::vector<graph::TreeHolder> cuts(const graph::Tree & A)
{
    std::vector<graph::TreeHolder> R;
    auto n=A.size();
    for(int i=0;i<n;i++) if(A.children(i).empty() || i==A.root && A.children(i).size()==1)
        R.push_back(cut(A,i));
    return R;
}

template<typename Mapper>
double expected_connectivity(const graph::TreeHolder &A,Mapper &M)
{
    auto n=A->size();
    if(n<=1)
        return n;
    if(M.count(A))
        return M[A];
    auto C=cuts(*A);
    double r=0;
    for(auto &B:C)
        r+=(expected_connectivity(B,M)+1)/n;
    return M[A]=r;
}

template<typename Mapper>
std::pair<double,size_t> expected_connectivity(const graph::TreeHolder &A)
{
    Mapper M;
    auto E=expected_connectivity(A,M);
    return std::make_pair(E,M.size());
}

int main()
{
    std::ios_base::sync_with_stdio(false);
    int n;
    std::cin >> n;
    graph::TreeHolder T1(n);
    for(int i=0;i<n-1;i++)
    {
        int u,v;
        std::cin >> u >> v;
        u--;
        v--;
        T1->connect(u,v);
    }
    T1->reRoot(0);
    T1->buildStatistics(graph::TreeStats::SIZE);
    auto t1=std::chrono::high_resolution_clock::now();
    auto [E,S]=expected_connectivity<TreeesHashMap>(T1);
    auto t2=std::chrono::high_resolution_clock::now();
    std::cout << E << std::endl;
    std::cout << "Size: " << S << std::endl;
    std::cout << "Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms" << std::endl;
}
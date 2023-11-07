#include "graph/tree/tree.h"
#include "graph/tree/isomorphism.h"
#include <iostream>

int main()
{
    std::ios_base::sync_with_stdio(false);
    int T;
    std::cin >> T;
    graph::IsoTreeEq treeEq;
    using TreesSet=std::set<graph::TreeHolder,graph::IsoTreeCmp>;
    for(int t=1;t<=T;t++)
    {
        TreesSet Z;
        int n;
        std::cin >> n;
        graph::TreeHolder T1(n), T2(n);
        for(int i=0;i<n-1;i++)
        {
            int u,v;
            std::cin >> u >> v;
            u--;
            v--;
            T1->connect(u,v);
        }
        for(int i=0;i<n-1;i++)
        {
            int u,v;
            std::cin >> u >> v;
            u--;
            v--;
            T2->connect(u,v);
        }
        T1->reRoot(0);
        T2->reRoot(0);
        T1->buildStatistics(graph::TreeStats::SIZE);
        T2->buildStatistics(graph::TreeStats::SIZE);
        Z.emplace(std::move(T1));
        Z.emplace(std::move(T2));
        auto result=(Z.size()==1?"YES":"NO");
        std::cout << result << '\n';
    }
}
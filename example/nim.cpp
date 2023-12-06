#include <iostream>
#include "graph/tree/isomorphism.h"

int main()
{
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int T;
    std::cin >> T;
    for(int t=1;t<=T;t++)
    {
        using Container=std::map<std::vector<int>,int>;
        std::unordered_set<cp::graph::TreeHolder,cp::graph::FastIsoTreeHash<Container>,cp::graph::FastIsoTreeEq<Container>> container;
        int n;
        std::cin >> n;
        cp::graph::TreeHolder T1(n),T2(n);
        for(auto T:{T1,T2}) for(int i=0;i<n-1;i++)
            {
                int u,v;
                std::cin >> u >> v;
                u--;
                v--;
                T->connect(u,v);
            }
        T1->reRoot(0);
        T2->reRoot(0);
        T1->buildStatistics(cp::graph::TreeStats::SIZE);
        T2->buildStatistics(cp::graph::TreeStats::SIZE);
        container.emplace(std::move(T1));
        container.emplace(std::move(T2));
        std::cout << (container.size()==1?"YES":"NO") << '\n';
    }
}
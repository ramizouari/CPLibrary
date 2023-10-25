//
// Created by ramizouari on 24/10/23.
//
#include <iostream>
#include "graph/tree/range_queries.h"

int main()
{
    int T;
    std::cin >> T;
    for(int t=1;t<=T;t++)
    {
        int n;
        std::cin >> n;
        HeavyLightTree<plus_t<int>,segment_tree<plus_t<int>>> H(n);
        std::vector<std::pair<int,int>> edges;
        for(int i=0;i<n-1;i++)
        {
            int u,v,w;
            std::cin >> u >> v >> w;
            u--,v--;
            H.connect(u,v,w);
            edges.emplace_back(u,v);
        }
        H.reRoot(0);
        H.buildStatistics();
        std::string query;
        std::cin >> query;
        while(query!="DONE")
        {
            if(query=="QUERY")
            {
                int u,v;
                std::cin >> u >> v;
                u--,v--;
                std::cout << H.query(u,v) << std::endl;
            }
            else if(query=="CHANGE")
            {
                int t,w;
                std::cin >> t >> w;
                t--;
                auto [u,v]=edges[t];
                H.update(u,v,w);
            }
            std::cin >> query;
        }
    }
}
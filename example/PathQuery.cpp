#include <iostream>
#include "graph/tree/range_queries.h"
#include <chrono>

int main()
{
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    int n,q;
    std::cin >> n >> q;
    CommutativeHeavyLightNodeTree<max_t<integer>,segment_tree<max_t<integer>>> H(n);
    std::vector<integer> A(n);
    for(auto &a:A)
        std::cin >> a;
    for(int i=0;i<n;i++)
        H.setWeight(i,A[i]);
    for(int i=0;i<n-1;i++)
    {
        int u,v;
        std::cin >> u >> v;
        u--;
        v--;
        H.connect(u,v);
    }
    H.build(0);
    auto t3=std::chrono::high_resolution_clock::now();
    for(int i=0;i<q;i++)
    {
        int b;
        std::cin >> b;
        if(b==1)
        {
            int s;
            integer x;
            std::cin >> s >> x;
            s--;
            H.update(s,x);
        }
        else
        {
            int u,v;
            std::cin >> u >> v;
            u--;
            v--;
            std::cout << H.query(u,v) << '\n';
        }
    }
    auto t4=std::chrono::high_resolution_clock::now();
    //std::cerr << "dt: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << std::endl;
}
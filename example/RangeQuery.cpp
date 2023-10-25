//
// Created by ramizouari on 25/10/23.
//
#include <iostream>
#include <data_structures/dynamic_range_queries.h>


int main()
{
    int n,q;
    std::cin >> n >> q;
    std::vector<integer> A(n);
    for(auto &a:A)
        std::cin >> a;
    auto F=std::make_shared<plus_t<integer>>();
    auto T=segment_tree<integer,void>(A,F);
    for(int i=0;i<q;i++)
    {
        int a,b;
        std::cin >> a >> b;
        std::cout << T.query(a,b) << '\n';
    }
}
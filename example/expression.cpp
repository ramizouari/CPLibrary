#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <numeric>

struct int_with_id
{
    int x,id;
    int_with_id(int _x,int _id):x(_x),id(_id){}
    int_with_id(int _x):x(_x){}
    std::weak_ordering operator<=>(const int_with_id &other) const
    {
        return x<=>other.x;
    }

    std::weak_ordering operator<=>(const int other) const
    {
        return x<=>other;
    }

    operator int() const
    {
        return id+1;
    }
};

using integer=long long;
int main()
{
    std::ios_base::sync_with_stdio(false);
    integer n,x;
    std::multiset<int_with_id> S;
    std::cin>>n >> x;
    for(integer i=0;i<n;i++)
    {
        integer x;
        std::cin>>x;
        S.emplace(x,i);
    }
    auto it1=S.begin();
    while(it1!=S.end())
    {
        auto it2=S.find(x-it1->x);
        if(it2!=S.end() && it1!=it2)
        {
            std::cout<<*it1<<" "<<*it2<<std::endl;
            return 0;
        }
        it1++;
    }
    std::cout<<"IMPOSSIBLE"<<std::endl;
}
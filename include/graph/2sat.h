//
// Created by ramizouari on 21/11/22.
//

#ifndef CPLIBRARY_2SAT_H
#define CPLIBRARY_2SAT_H

#include "graph.h"
#include <stack>
#include <ostream>
#include <optional>

constexpr int Id(int x)
{
    return x<<1;
}

constexpr int Not(int x)
{
    return (x<<1)|1;
}

constexpr int Var(int x,bool r)
{
    return r?Id(x):Not(x);
}

constexpr int BaseVariable(int r)
{
    return r>>1;
}

constexpr bool isNegation(int r)
{
    return r&1;
}

constexpr bool isAffirmative(int r)
{
    return !isNegation(r);
}


class SatisfiabilityGraph:public Graph
{
    std::vector<int> assumptions;
    using Graph::n;
    using Graph::adjacencyList;
    using Graph::reverseList;
    int d;
public:
    explicit SatisfiabilityGraph(int d):Graph(2*d),d(d){}
    void connect(int a,int b,bool r1,bool r2)
    {
        adjacencyList[Var(a,r1)].push_back(Var(b,r2));
        reverseList[Var(b,r2)].push_back(Var(a,r1));
    }

    void addClause(int a,int b,bool r1, bool r2)
    {
        connect(a,b,r1,r2);
        connect(b,a,!r2,!r1);
    }

    void addOr(int a,int b)
    {
        addClause(a,b,false,true);
    }

    void addNand(int a,int b)
    {
        addClause(a,b,true,false);
    }

    void addImplication(int a,int b)
    {
        addClause(a,b,true,true);
    }

    void addEquivalence(int a,int b)
    {
        addImplication(a,b);
        addImplication(b,a);
    }

    void addExclusion(int a,int b)
    {
        addOr(a,b);
        addNand(a,b);
    }


    void addAssumption(int a,bool r)
    {
        addClause(a,a,!r,r);
        assumptions.push_back(Var(a,r));
    }

    bool satisfiable()
    {
        return satisfy().has_value();
    }

    std::optional<std::vector<bool>> satisfy()
    {
        std::stack<int> Q;
        std::vector<bool> visited(2*d);
        std::vector<bool> solution(d),assigned(d);
        auto &&[components,componentId,classes,topologicalOrder]=getConnectedComponentsWithMetaData();
        std::vector<bool> componentAssigned(components.size());

        std::vector<bool> visitedComponents(2*d);
        for(auto &component:components)
        {
            for(auto &u:component)
                visitedComponents[u]=true;
            for(auto &u:component)
            {
                int i= BaseVariable(u);
                if(visitedComponents[Id(i)] && visitedComponents[Not(i)])
                    return std::nullopt;
            }
            for(auto &u:component)
                visitedComponents[u]=false;
        }
        for(int i=topologicalOrder.size()-1;i>=0;i--) if(!componentAssigned[componentId[topologicalOrder[i]]])
            {
                componentAssigned[i]=true;
                auto &component=components[componentId[topologicalOrder[i]]];
                bool unassigned = std::none_of(component.begin(),component.end(),[&assigned](const auto &s){return assigned[BaseVariable(s)];});
                if(unassigned) for(auto c:component)
                    {
                        assigned[BaseVariable(c)]=true;
                        solution[BaseVariable(c)]=isAffirmative(c);
                    }
            }
        return solution;
    }

    void print(std::ostream &H,int i) const
    {
        if(i&1)
            H << "¬";
        H << "X[" << BaseVariable(i)  << "]";
    }

    void print(std::ostream &H) const
    {
        for(auto u:assumptions) {
            print(H, u);
            H << '\n';
        }
        for(int i=0;i<2*d;i++)
        {
            for(auto u:adjacencyList[i]) {
                print(H,i);
                H << " -> ";
                print(H,u);
                H << '\n';
            }
        }
    }
};

#endif //CPLIBRARY_2SAT_H

//
// Created by ramizouari on 21/11/22.
//

#ifndef CPLIBRARY_2SAT_H
#define CPLIBRARY_2SAT_H

#include "graph.h"
#include <stack>
#include <ostream>
#include <optional>
namespace cp::sat
{
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


    struct Sat2 : graph::Graph
    {
        int d;
        using Graph=cp::graph::Graph;
        explicit Sat2(int d) : Graph(2*d),d(d) {}

        // Get a solution, or an empty vector if not satisfiable
        std::optional<std::vector<bool>> satisfy()
        {
            std::vector<bool> vis(2 * d), sol(d), assigned(d);
            auto &&[components,cmpId, _] = getSCComponents();
            for (int i = 0; i < n; i += 2) {
                if (cmpId[i] == cmpId[i^1])
                    return std::nullopt;
                sol[i / 2] = cmpId[i] > cmpId[i^1];
            }
            return sol;
        }

        // Connect two nodes in a graph
        void connect(int a, int b, bool r1, bool r2)
        {
            Graph::connect(2*a^r1,2*b^r2);
        }

        // Add a 2-SAT clause.
        void addClause(int a, int b, bool r1, bool r2) {
            connect(a, b, !r1, r2);
            connect(b, a, !r2, r1);
        }

        // a ∨ b
        void addOr(int a, int b) {
            addClause(a, b, false, false);
        }

        // ⅂a ∨ ⅂b
        void addNand(int a, int b) {
            addClause(a, b, true, true);
        }

        // a => b. Equivalent to: ⅂a ∨ b
        void addImplication(int a, int b) {
            addClause(a, b, true, false);
        }

        // a <=> b. Equivalent to a => b and b => a
        void addEquivalence(int a, int b) {
            addImplication(a, b);
            addImplication(b, a);
        }

        // a ∨ b, but not a ∧ b
        void addExclusion(int a, int b) {
            addOr(a, b);
            addNand(a, b);
        }

        // Assume that a has value r.
        void addAssumption(int a, bool r) {
            addClause(a, a, !r, r);
        }

        bool satisfiable()
        {
            return !satisfy();
        }
    };
}

#endif //CPLIBRARY_2SAT_H

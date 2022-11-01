//
// Created by ramizouari on 31/10/22.
//

#include <unordered_set>
#include <queue>
#include "formal_verification/abstract_printer.h"

namespace FormalSpecification::Printer {
        abstract_printer::abstract_printer(std::ostream &S):ostream(S) {

        }

    void mermaid_printer::print(const dynamic_decision_tree &T)
    {
        auto root=T.root;
        auto &H=ostream;
        std::unordered_set<FormalSpecification::node_t*> S;
        std::unordered_map<FormalSpecification::node_t*,std::uint64_t> id;
        std::unordered_map<std::uint64_t ,std::uint64_t> C;
        std::queue<FormalSpecification::node_t*> Q;
        Q.push(root);
        S.emplace(root);
        id[root]=0;
        C[root->symbol]=1;
        H << "flowchart TB" << '\n';
        if(root->is_leaf())
            H << "\t" << (root->symbol==Truth::True?"True":"False") << '\n';
        else while(!Q.empty())
        {
            auto u=Q.front();
            Q.pop();
            if(u->is_leaf())
                continue;

            if(!S.contains(u->left)) {
                Q.push(u->left);
                S.emplace(u->left);
                auto r=C[u->left->symbol]+1;
                id[u->left]=r;
                C[u->left->symbol]=r;
            }
            if(!S.contains(u->right)) {
                Q.push(u->right);
                S.emplace(u->right);
                auto r=C[u->right->symbol]+1;
                id[u->right]=r;
                C[u->right->symbol]=r;
            }

            H << "\t" << u->symbol << 'I' << id[u] << '[' << u->symbol <<  ']' << " --False--> ";
            if(!u->left->is_leaf())
                H << u->left->symbol << 'I' << id[u->left];
            else
                H << (u->left->symbol==FormalSpecification::True?"True":"False");
            H << '\n';
            H << "\t" << u->symbol << 'I' << id[u] << " --True--> ";
            if(!u->right->is_leaf())
                H << u->right->symbol << 'I' << id[u->right];
            else
                H << (u->right->symbol==FormalSpecification::True?"True":"False");
            H << '\n';

        }
    }

    abstract_printer& operator<<(abstract_printer &P, const dynamic_decision_tree &T) {
            P.print(T);
            return P;
    }

    void graphml_printer::print(const dynamic_decision_tree &T)
    {
        auto root=T.root;
        auto &H=ostream;
        std::unordered_set<FormalSpecification::node_t*> S;
        std::unordered_map<FormalSpecification::node_t*,std::uint64_t> id;
        std::unordered_map<std::uint64_t ,std::uint64_t> C;
        std::queue<FormalSpecification::node_t*> Q;
        Q.push(root);
        S.emplace(root);
        id[root]=0;
        C[root->symbol]=1;
        H << "flowchart TB" << '\n';
        while(!Q.empty())
        {
            auto u=Q.front();
            Q.pop();
            if(u->is_leaf())
                continue;

            if(!S.contains(u->left)) {
                Q.push(u->left);
                S.emplace(u->left);
                auto r=C[u->left->symbol]+1;
                id[u->left]=r;
                C[u->left->symbol]=r;
            }
            if(!S.contains(u->right)) {
                Q.push(u->right);
                S.emplace(u->right);
                auto r=C[u->right->symbol]+1;
                id[u->right]=r;
                C[u->right->symbol]=r;
            }

            H << "\t" << u->symbol << 'I' << id[u] << '[' << u->symbol <<  ']' << " --False--> ";
            if(!u->left->is_leaf())
                H << u->left->symbol << 'I' << id[u->left];
            else
                H << (u->left->symbol==FormalSpecification::True?"True":"False");
            H << '\n';
            H << "\t" << u->symbol << 'I' << id[u] << " --True--> ";
            if(!u->right->is_leaf())
                H << u->right->symbol << 'I' << id[u->right];
            else
                H << (u->right->symbol==FormalSpecification::True?"True":"False");
            H << '\n';

        }
    }
} // Printer
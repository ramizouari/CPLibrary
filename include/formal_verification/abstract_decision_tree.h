//
// Created by ramizouari on 29/10/22.
//

#ifndef VFSF_ABSTRACT_DECISION_TREE_H
#define VFSF_ABSTRACT_DECISION_TREE_H

#include "node.h"

namespace FormalSpecification
{
    struct abstract_decision_tree
    {
        node_t *root;
        int variables;
        abstract_decision_tree(int _variables,node_t* _root=nullptr):root(_root),variables(_variables){}
        virtual ~abstract_decision_tree()=0;
        std::uint64_t count_truths() const;
    };

}

#endif //VFSF_ABSTRACT_DECISION_TREE_H

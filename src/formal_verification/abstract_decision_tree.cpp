//
// Created by ramizouari on 29/10/22.
//


#include "formal_verification/abstract_decision_tree.h"

namespace FormalSpecification
{
    abstract_decision_tree::~abstract_decision_tree(){}

    std::uint64_t abstract_decision_tree::count_truths() const {
        return FormalSpecification::count_truths(root,variables);
    }
}
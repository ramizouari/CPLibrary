//
// Created by ramizouari on 27/05/22.
//
#include "parser/LRParserBuilder.h"

namespace cp::parser
{
    std::unordered_set<LR0Item> LR0ParserBuilder::closure(const std::unordered_set<LR0Item> &items) {
        std::unordered_set<LR0Item> newItems;
        std::stack<std::uint64_t> stack;
        for(const auto& item:items)
        {
            newItems.insert(item);
            if(item.dot < rules[item.ruleId].right.size())
            {
                if(isTerminal[rules[item.ruleId].right[item.dot].id])
                    continue;
                stack.push(rules[item.ruleId].right[item.dot].id);
            }
            while(!stack.empty())
            {
                auto id=stack.top();
                stack.pop();
                for(const auto & ruleId : ruleIdsBySymbolId[id])
                {
                    auto [_,inserted]=newItems.insert(LR0Item{ruleId,0});
                    if(inserted && !rules[ruleId].right.empty() && !isTerminal[rules[ruleId].right.front().id])
                        stack.push(rules[ruleId].right.front().id);
                }
            }
        }
        return newItems;
    }

    LR0ParserBuilder &LR0ParserBuilder::build() &
    {
        Symbol augmentedSymbol=addSymbol("");
        Grammar::addRule(Rule(augmentedSymbol,{axiom}));
        setAxiom(augmentedSymbol);
        auto augmentedRuleId=rules.size()-1;
        ruleIdsBySymbolId.resize(symbols.size());
        isTerminal.resize(symbols.size(),true);
        for(const auto & rule : rules)
            isTerminal[rule.left.id]=false;
        LR0Item I0{rules.size()-1,0};
        lr0Items.push_back(std::move(closure({I0})));
        lr0ItemsIds.emplace(lr0Items.back(),lr0Items.size()-1);
        int k=0,n=1;
        do
        {
            std::unordered_map<std::uint64_t,std::unordered_set<LR0Item>> nextItemsSets;
            for(const auto & item : lr0Items[k])
            {
                if(item.dot<rules[item.ruleId].right.size())
                    nextItemsSets[rules[item.ruleId].right[item.dot].id].emplace(item.ruleId,item.dot+1);
            }
            for(auto &[symbolId,itemsSet]: nextItemsSets)
            {
                itemsSet = std::move(closure(itemsSet));
                auto it=lr0ItemsIds.find(itemsSet);
                if(it==lr0ItemsIds.end())
                {
                    lr0Items.push_back(std::move(itemsSet));
                    std::tie(it,std::ignore)=lr0ItemsIds.emplace(lr0Items.back(),lr0Items.size()-1);
                    n++;
                }
                gotoIds.emplace(std::make_pair(k,symbolId),Action(Action::Shift,it->second));
            }
        } while(++k<n);
        for(int i=0;i<n;i++) for(const auto &item: lr0Items[i]) if(rules[item.ruleId].right.size()==item.dot)
        {
            for (auto &F: symbols) if (isTerminal[F.id])
            {
                if (item.ruleId == augmentedRuleId)
                    gotoIds.emplace(std::make_pair(i, SpecialCharacter::EndOfString),
                                    Action(Action::Accept, 0));
                else gotoIds.emplace(std::make_pair(i, F.id), Action(Action::Reduce, item.ruleId));
            }
            gotoIds.emplace(std::make_pair(i, SpecialCharacter::EndOfString), Action(Action::Reduce, item.ruleId));
        }
        return *this;
    }

    StatefulLR0ParserBuilder::StatefulLR0ParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds)
    {

    }

}
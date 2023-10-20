//
// Created by ramizouari on 27/05/22.
//
#include <queue>
#include "parser/LRParserBuilder.h"

namespace parser
{
    std::unordered_set<LALR1Item> LALRParserBuilder::closure(const std::unordered_set<LALR1Item> &items)
    {
        std::unordered_set<LALR1Item> newItems;
        struct context
        {
            std::uint64_t symbolId;
            std::unordered_set<std::uint64_t> firstIds;
        };
        std::stack<context> stack;
        for(const auto& item:items)
        {
            auto [it,inserted] = newItems.insert(item);
            if(!inserted)
                it->lookahead.insert(item.lookahead.begin(),item.lookahead.end());

            if(item.dot < rules[item.ruleId].right.size())
            {
                if(isTerminal[rules[item.ruleId].right[item.dot].id])
                    continue;
                context ctx;
                ctx.symbolId=rules[item.ruleId].right[item.dot].id;
                if(item.dot+1==rules[item.ruleId].right.size())
                    ctx.firstIds.insert(item.lookahead.begin(),item.lookahead.end());
                else for(auto k=item.dot+1;k<rules[item.ruleId].right.size();k++)
                {
                    ctx.firstIds.insert(firstIds[rules[item.ruleId].right[k].id].begin(),firstIds[rules[item.ruleId].right[k].id].end());
                    if(!firstIds[rules[item.ruleId].right[k].id].contains(SpecialCharacter::Epsilon))
                        break;
                }
                ctx.firstIds.erase(SpecialCharacter::Epsilon);
                stack.push(ctx);
            }
            while(!stack.empty())
            {
                auto [id,lookaheads]=stack.top();
                stack.pop();
                bool anyInsertion=false;
                for(const auto & ruleId : ruleIdsBySymbolId[id])
                {
                    context ctx;
                    if(!rules[ruleId].right.empty())
                        ctx.symbolId = rules[ruleId].right.front().id;
                    auto itemIterator=newItems.find(LALR1Item(ruleId, 0, {}));
                    if(itemIterator==newItems.end())
                    {
                        newItems.insert(LALR1Item(ruleId, 0, lookaheads));
                        anyInsertion=true;
                    }
                    else for(auto lookahead:lookaheads) if(!itemIterator->lookahead.contains(lookahead))
                    {
                        itemIterator->lookahead.insert(lookahead);
                        anyInsertion=true;
                    }
                    if (anyInsertion && !rules[ruleId].right.empty() && !isTerminal[rules[ruleId].right.front().id])
                    {
                        int k;
                        for (k = 1; k < rules[ruleId].right.size(); k++)
                        {
                            ctx.firstIds.insert(firstIds[rules[ruleId].right[k].id].begin(),
                                                firstIds[rules[ruleId].right[k].id].end());
                            if (!firstIds[rules[ruleId].right[k].id].contains(SpecialCharacter::Epsilon))
                                break;
                        }
                        if(k==rules[ruleId].right.size())
                            ctx.firstIds.insert(lookaheads.begin(),lookaheads.end());
                        ctx.firstIds.erase(SpecialCharacter::Epsilon);
                    }
                    if(anyInsertion)
                        stack.push(ctx);
                }
            }
        }
        return newItems;
    }

    LALRParserBuilder &LALRParserBuilder::build() &{
        Symbol augmentedSymbol=addSymbol("");
        Grammar::addRule(Rule(augmentedSymbol,{axiom}));
        setAxiom(augmentedSymbol);
        auto augmentedRuleId=rules.size()-1;
        ruleIdsBySymbolId.resize(symbols.size());
        buildFirst();
        buildFollow();
        isTerminal.resize(symbols.size(),true);

        for(const auto & rule : rules)
            isTerminal[rule.left.id]=false;
        LALR1Item I0{rules.size()-1,0,{SpecialCharacter::EndOfString}};
        lr1Items.push_back(std::move(closure({I0})));
        lr1ItemsIds.emplace(lr1Items.back(), lr1Items.size() - 1);
        /*
         * Queue of Item-sets that need to be updated.
         * */
        std::queue<int> updates;
        updates.emplace(0);
        int n=1;
        do
        {
            auto k = updates.front();
            updates.pop();
            std::unordered_map<std::uint64_t,std::unordered_set<LALR1Item>> nextItemsSets;
            auto update=lr1ItemsIds.find(lr1Items[k]);
            lr1Items[k] = update->first;
            for(const auto & item : lr1Items[k])
            {
                if(item.dot<rules[item.ruleId].right.size())
                {
                    auto [it,inserted] = nextItemsSets[rules[item.ruleId].right[item.dot].id].emplace(item.ruleId, item.dot + 1, item.lookahead);
                    if(!inserted)
                        throw std::runtime_error("Debugging: duplicate item in LALR1 closure");
                }
            }
            for(auto &[symbolId,itemsSet]: nextItemsSets)
            {
                itemsSet = std::move(closure(itemsSet));
                auto it=lr1ItemsIds.find(itemsSet);
                if(it == lr1ItemsIds.end())
                {
                    lr1Items.push_back(std::move(itemsSet));
                    std::tie(it,std::ignore)=lr1ItemsIds.emplace(lr1Items.back(), lr1Items.size() - 1);
                    updates.emplace(it->second);
                    n++;
                }
                //Merge lookaheads
                else
                {
                    bool insertion=false;
                    for (auto &item: it->first)
                    {
                        auto item2 = itemsSet.find(item);
                        if (item2 != itemsSet.end()) for (auto lookahead : item2->lookahead) if (!item.lookahead.contains(lookahead))
                        {
                            item.lookahead.insert(lookahead);
                            insertion=true;
                        }
                    }
                    if(insertion)
                        updates.emplace(it->second);
                }
                gotoIds.emplace(std::make_pair(k,symbolId),Action(Action::Shift,it->second));
            }
        } while(!updates.empty());
        for(int i=0;i<n;i++) for(const auto &item: lr1Items[i]) if(rules[item.ruleId].right.size() == item.dot)
        {
            if (item.ruleId == augmentedRuleId)
                gotoIds.emplace(std::make_pair(i, SpecialCharacter::EndOfString), Action(Action::Accept, 0));
            else for(auto lookahead : item.lookahead)
                    gotoIds.emplace(std::make_pair(i, lookahead), Action(Action::Reduce, item.ruleId));
        }
        return *this;
    }


    StatefulLALRParserBuilder::StatefulLALRParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds)
    {

    }
}
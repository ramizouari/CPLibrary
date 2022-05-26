//
// Created by ramizouari on 26/05/22.
//

#include "parser/LRParser.h"

namespace parser {


    SLRParser& SLRParser::build() &
    {
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
            for(auto &F:followIds[rules[item.ruleId].left.id])
            {
                if (item.ruleId == augmentedRuleId)
                    gotoIds.emplace(std::make_pair(i,SpecialCharacter::EndOfString),Action(Action::Accept,0));
                else    gotoIds.emplace(std::make_pair(i,F),Action(Action::Reduce,item.ruleId));
            }
        return *this;
    }

    bool SLRParser::parse(const std::string &S)
    {
        std::stack<std::pair<std::uint64_t,std::uint64_t>> stack;
        stack.emplace(0,-1);
        char character[2]{0,0};
        for(int i=0;i<S.size();)
        {
            auto &s=S[i];
            auto [state,ruleId]=stack.top();
            character[0]=s;
            auto symbolId=symbolMap.at(character);
            auto it=gotoIds.find(std::make_pair(state,symbolId));
            if(it==gotoIds.end())
                return false;
            auto [type,value]=it->second;
            switch(type)
            {
                case Action::Shift:
                    i++;
                    stack.emplace(value,symbolId);
                    break;
                case Action::Reduce:
                {
                    if (stack.size() < rules[value].right.size())
                        return false;
                    for (int j = 0; j < rules[value].right.size(); j++)
                        stack.pop();
                    auto newState = gotoIds.find(std::make_pair(stack.top().first, rules[value].left.id));
                    if (newState == gotoIds.end())
                        return false;
                    stack.emplace(newState->second.value,rules[value].left.id);
                    break;
                }
                case Action::Accept:
                    return true;

            }
        }
        while(!stack.empty())
        {
            auto [state,ruleId]=stack.top();
            auto it=gotoIds.find(std::make_pair(state,SpecialCharacter::EndOfString));
            if(it==gotoIds.end())
                return false;
            auto [type,value]=it->second;
            switch(type)
            {
                case Action::Reduce:
                {
                    if (stack.size() < rules[value].right.size())
                        return false;
                    for (int j = 0; j < rules[value].right.size(); j++)
                        stack.pop();
                    auto newState = gotoIds.find(std::make_pair(stack.top().first, rules[value].left.id));
                    if (newState == gotoIds.end())
                        return false;
                    stack.emplace(newState->second.value,rules[value].left.id);
                    break;
                }
                case Action::Accept:
                    return true;
            }
        }
        return false;
    }

    std::unordered_set<LR0Item> SLRParser::closure(const std::unordered_set<LR0Item> &items)
    {
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


    LR0Item::LR0Item(std::uint64_t ruleId, std::uint64_t dot) : ruleId(ruleId), dot(dot)
    {
        initHash();
    }

    std::uint64_t LR0Item::initHash() {
        std::uint64_t hash=ruleId;
        hash^=hasher(dot);
        return hash;
    }

    LRFamily::Action::Action(LRFamily::Action::Type type, std::uint64_t value) : type(type), value(value)
    {
    }


    LRParser &LRParser::build() &
    {
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
        LR1Item I0{rules.size()-1,0,{SpecialCharacter::EndOfString}};
        lr1Items.push_back(std::move(closure({I0})));
        lr1ItemsIds.emplace(lr1Items.back(), lr1Items.size() - 1);
        int k=0,n=1;
        do
        {
            std::unordered_map<std::uint64_t,std::unordered_set<LR1Item>> nextItemsSets;
            for(const auto & item : lr1Items[k])
            {
                if(item.dot<rules[item.ruleId].right.size())
                    nextItemsSets[rules[item.ruleId].right[item.dot].id].emplace(item.ruleId,item.dot+1,item.lookahead);
            }
            for(auto &[symbolId,itemsSet]: nextItemsSets)
            {
                itemsSet = std::move(closure(itemsSet));
                auto it=lr1ItemsIds.find(itemsSet);
                if(it == lr1ItemsIds.end())
                {
                    lr1Items.push_back(std::move(itemsSet));
                    std::tie(it,std::ignore)=lr1ItemsIds.emplace(lr1Items.back(), lr1Items.size() - 1);
                    n++;
                }
                gotoIds.emplace(std::make_pair(k,symbolId),Action(Action::Shift,it->second));
            }
        } while(++k<n);
        for(int i=0;i<n;i++) for(const auto &item: lr1Items[i]) if(rules[item.ruleId].right.size() == item.dot)
                {
                    if (item.ruleId == augmentedRuleId)
                        gotoIds.emplace(std::make_pair(i, SpecialCharacter::EndOfString), Action(Action::Accept, 0));
                    else gotoIds.emplace(std::make_pair(i, item.lookahead), Action(Action::Reduce, item.ruleId));
                }
        return *this;
    }

    std::unordered_set<LR1Item> LRParser::closure(const std::unordered_set<LR1Item> &items) {
        std::unordered_set<LR1Item> newItems;
        struct context
        {
            std::uint64_t symbolId;
            std::unordered_set<std::uint64_t> firstIds;
        };
        std::stack<context> stack;
        for(const auto& item:items)
        {
            newItems.insert(item);
            if(item.dot < rules[item.ruleId].right.size())
            {
                if(isTerminal[rules[item.ruleId].right[item.dot].id])
                    continue;
                context ctx;
                ctx.symbolId=rules[item.ruleId].right[item.dot].id;
                if(item.dot+1==rules[item.ruleId].right.size())
                    ctx.firstIds.insert(item.lookahead);
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
                    for (auto lookahead: lookaheads)
                    {
                        auto [newItemIterator, inserted] = newItems.insert(LR1Item{ruleId, 0, lookahead});
                        anyInsertion = anyInsertion || inserted;
                        if (inserted && !rules[ruleId].right.empty() && !isTerminal[rules[ruleId].right.front().id])
                        {
                            if (rules[ruleId].right.size() == 1)
                                ctx.firstIds.insert(lookahead);
                            else
                                for (auto k = 1; k < rules[ruleId].right.size(); k++)
                                {
                                    ctx.firstIds.insert(firstIds[rules[ruleId].right[k].id].begin(),
                                                        firstIds[rules[ruleId].right[k].id].end());
                                    if (!firstIds[rules[ruleId].right[k].id].contains(SpecialCharacter::Epsilon))
                                        break;
                                }
                            ctx.firstIds.erase(SpecialCharacter::Epsilon);
                        }
                    }
                    if(anyInsertion)
                        stack.push(ctx);
                }
            }
        }
        return newItems;
    }

    bool LRParser::parse(const std::string &S) {
        std::stack<std::pair<std::uint64_t,std::uint64_t>> stack;
        stack.emplace(0,-1);
        char character[2]{0,0};
        for(int i=0;i<S.size();)
        {
            auto &s=S[i];
            auto [state,ruleId]=stack.top();
            character[0]=s;
            auto symbolId=symbolMap.at(character);
            auto it=gotoIds.find(std::make_pair(state,symbolId));
            if(it==gotoIds.end())
                return false;
            auto [type,value]=it->second;
            switch(type)
            {
                case Action::Shift:
                    i++;
                    stack.emplace(value,symbolId);
                    break;
                case Action::Reduce:
                {
                    if (stack.size() < rules[value].right.size())
                        return false;
                    for (int j = 0; j < rules[value].right.size(); j++)
                        stack.pop();
                    auto newState = gotoIds.find(std::make_pair(stack.top().first, rules[value].left.id));
                    if (newState == gotoIds.end())
                        return false;
                    stack.emplace(newState->second.value,rules[value].left.id);
                    break;
                }
                case Action::Accept:
                    return true;

            }
        }
        while(!stack.empty())
        {
            auto [state,ruleId]=stack.top();
            auto it=gotoIds.find(std::make_pair(state,SpecialCharacter::EndOfString));
            if(it==gotoIds.end())
                return false;
            auto [type,value]=it->second;
            switch(type)
            {
                case Action::Reduce:
                {
                    if (stack.size() < rules[value].right.size())
                        return false;
                    for (int j = 0; j < rules[value].right.size(); j++)
                        stack.pop();
                    auto newState = gotoIds.find(std::make_pair(stack.top().first, rules[value].left.id));
                    if (newState == gotoIds.end())
                        return false;
                    stack.emplace(newState->second.value,rules[value].left.id);
                    break;
                }
                case Action::Accept:
                    return true;
            }
        }
        return false;
    }

    LR1Item::LR1Item(std::uint64_t ruleId, std::uint64_t dot, std::uint64_t lookahead) : ruleId(ruleId), dot(dot), lookahead(lookahead)
    {
    }

} // parser
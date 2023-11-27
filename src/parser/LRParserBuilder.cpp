//
// Created by ramizouari on 26/05/22.
//

#include <iomanip>
#include "parser/LRParserBuilder.h"
#include "parser/StatefulParser.h"

namespace cp::parser {

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


    LRParserBuilder &LRParserBuilder::build() &
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
            [[unlikely]]
            if (item.ruleId == augmentedRuleId)
                gotoIds.emplace(std::make_pair(i, SpecialCharacter::EndOfString), Action(Action::Accept, 0));
            else gotoIds.emplace(std::make_pair(i, item.lookahead), Action(Action::Reduce, item.ruleId));
        }
        return *this;
    }

    std::unordered_set<LR1Item> LRParserBuilder::closure(const std::unordered_set<LR1Item> &items) {
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
                           int k;
                            for (k = 1; k < rules[ruleId].right.size(); k++)
                            {
                                ctx.firstIds.insert(firstIds[rules[ruleId].right[k].id].begin(),
                                                    firstIds[rules[ruleId].right[k].id].end());
                                if (!firstIds[rules[ruleId].right[k].id].contains(SpecialCharacter::Epsilon))
                                    break;
                            }
                            if(k==rules[ruleId].right.size())
                                ctx.firstIds.insert(lookahead);
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

    bool ShiftReduceParser::parse(const std::string &S) {
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

    LRFamily::LALR1ItemContext::LALR1ItemContext(std::uint64_t id, std::unordered_set<std::uint64_t> &&lookaheads):id(id),lookaheads(std::move(lookaheads))
    {
    }

    LRFamily::LALR1ItemContext::LALR1ItemContext(std::uint64_t id,
                                                 const std::unordered_set<std::uint64_t> &lookaheads) : id(id),
                                                                                                       lookaheads(lookaheads)
                                                                                                       {

                                                                                                       }
    LALR1Item::LALR1Item(std::uint64_t ruleId, std::uint64_t dot,const std::unordered_set<std::uint64_t> &lookahead) : ruleId(ruleId), dot(dot), lookahead(lookahead)
    {
    }

    LALR1Item::LALR1Item(std::uint64_t ruleId, std::uint64_t dot, std::unordered_set<std::uint64_t> &&lookahead) : ruleId(ruleId), dot(dot), lookahead(std::move(lookahead))
    {
    }

    bool LALR1Item::operator==(const LALR1Item &other) const {
        return ruleId == other.ruleId && dot == other.dot;
    }


    void ShiftReduceParser::printTable(std::ostream &H) const{
        std::vector<std::vector<Action>> table;
        auto mapSymbol=[count=symbolId.size()-1](auto x) -> decltype(x){
            if(x==SpecialCharacter::EndOfString)
                return count;
            return x;
        };
        for(auto [cell,action]:gotoIds)
        {
            if(table.size()<cell.first+1)
                table.resize(cell.first+1,std::vector<Action>(symbolMap.size(),Action{Action::Error,0}));
           table[cell.first][mapSymbol(cell.second)]=action;
        }
        std::vector<std::string> symbolNames(symbolMap.size());
        for(auto &[symbol,id]:symbolMap)
        {
            symbolNames[id] = symbol;
            if(symbol == "\n")
                symbolNames[id] = "\\n";
        }
        symbolNames.back()="$";
        H << std::setw(5) << "" <<"\t";
        for(const auto& name:symbolNames)
            H << std::setw(5) << name << "\t";
        H << std::endl;
        for(int state=0;state<table.size();state++)
        {
            H  << std::setw(5) << "S" << state << '\t';
            for (int symbol = 0; symbol < table[state].size(); symbol++)
            {
                switch (table[state][mapSymbol(symbol)].type)
                {
                    case Action::Shift:
                        H << std::setw(5) << "s" + std::to_string(table[state][mapSymbol(symbol)].value) << '\t';
                        break;
                    case Action::Reduce:
                        H << std::setw(5) << "r" << std::to_string(table[state][mapSymbol(symbol)].value) << '\t';
                        break;
                    case Action::Accept:
                        H << std::setw(5) << "acc" << '\t';
                        break;
                    case Action::Error:
                    H << std::setw(5) << "err" << '\t';
                }
            }
            H << std::endl;
        }
    }



    std::shared_ptr<Variable> StatefulShiftReduceParser::evaluate(const std::string &S)
    {
        std::stack<std::tuple<std::uint64_t,std::uint64_t,std::shared_ptr<Variable>>> stack;
        stack.emplace(0,-1, nullptr);
        char character[2]{0,0};
        for(int i=0;i<S.size();)
        {
            auto &s=S[i];
            auto [state,ruleId,variable]=stack.top();
            character[0]=s;
            auto symbolId=symbolMap.at(character);
            auto it=gotoIds.find(std::make_pair(state,symbolId));
            if(it==gotoIds.end())
                return nullptr;
            auto [type,value]=it->second;
            switch(type)
            {
                case Action::Shift:
                    i++;
                    stack.emplace(value,symbolId,variable);
                    break;
                case Action::Reduce:
                {
                    std::vector<std::shared_ptr<Variable>> reduced_variables(rules[value].right.size());
                    if (stack.size() < rules[value].right.size())
                        return nullptr;
                    for (int j = 0; j < rules[value].right.size(); j++)
                    {
                        auto [_1,_2,var] = stack.top();
                        stack.pop();
                        reduced_variables[rules[value].right.size()-j-1]=var;
                    }
                    auto newState = gotoIds.find(std::make_pair(std::get<0>(stack.top()), rules[value].left.id));
                    if (newState == gotoIds.end())
                        return nullptr;
                    stack.emplace(newState->second.value,rules[value].left.id,reducers[value]->combine(reduced_variables));
                    break;
                }
                case Action::Accept:
                    return variable;
            }
        }
        while(!stack.empty())
        {
            auto [state,ruleId,variable]=stack.top();
            auto it=gotoIds.find(std::make_pair(state,SpecialCharacter::EndOfString));
            if(it==gotoIds.end())
                return nullptr;
            auto [type,value]=it->second;
            switch(type)
            {
                case Action::Reduce:
                {
                    if (stack.size() < rules[value].right.size())
                        return nullptr;
                    std::vector<std::shared_ptr<Variable>> reduced_variables(rules[value].right.size());
                    for (int j = 0; j < rules[value].right.size(); j++)
                    {
                        auto [_1,_2,var] = stack.top();
                        stack.pop();
                        reduced_variables[rules[value].right.size()-j-1]=var;
                    }
                    auto newState = gotoIds.find(std::make_pair(std::get<0>(stack.top()), rules[value].left.id));
                    if (newState == gotoIds.end())
                        return nullptr;
                    stack.emplace(newState->second.value,rules[value].left.id,reducers[value]->combine(reduced_variables));
                    break;
                }
                case Action::Accept:
                    return variable;
            }
        }
        return nullptr;
    }

    StatefulShiftReduceParser::StatefulShiftReduceParser(StatefulShiftReduceParser::IdsMapType &gotoIds) : gotoIds(gotoIds) {

    }

    StatefulLRParserBuilder::StatefulLRParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds){

    }

} // parser
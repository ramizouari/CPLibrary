//
// Created by ramizouari on 20/10/23.
//
#include "parser/LRLogic.h"

using namespace parser;

bool LRLogic::parse(const std::string &S)
{
    std::stack<std::tuple<std::uint64_t,std::uint64_t,std::shared_ptr<Variable>>> stack;
    stack.emplace(0,-1, nullptr);
    char character[2]{0,0};
    result = nullptr;
    for(int i=0;i<S.size();)
    {
        auto &s=S[i];
        auto [state,ruleId,variable]=stack.top();
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
                stack.emplace(value,symbolId,variable);
                break;
            case Action::Reduce:
            {
                std::vector<std::shared_ptr<Variable>> reduced_variables(rules[value].right.size());
                if (stack.size() < rules[value].right.size())
                    return false;
                for (int j = 0; j < rules[value].right.size(); j++)
                {
                    auto [_1,_2,var] = stack.top();
                    stack.pop();
                    reduced_variables[rules[value].right.size()-j-1]=var;
                }
                auto newState = gotoIds.find(std::make_pair(std::get<0>(stack.top()), rules[value].left.id));
                if (newState == gotoIds.end())
                    return false;
                stack.emplace(newState->second.value,rules[value].left.id,reducers[value]->combine(reduced_variables));
                break;
            }
            case Action::Accept:
            {
                result=variable;
                return true;
            }

        }
    }
    while(!stack.empty())
    {
        auto [state,ruleId,variable]=stack.top();
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
                std::vector<std::shared_ptr<Variable>> reduced_variables(rules[value].right.size());
                for (int j = 0; j < rules[value].right.size(); j++)
                {
                    auto [_1,_2,var] = stack.top();
                    stack.pop();
                    reduced_variables[rules[value].right.size()-j-1]=var;
                }
                auto newState = gotoIds.find(std::make_pair(std::get<0>(stack.top()), rules[value].left.id));
                if (newState == gotoIds.end())
                    return false;
                stack.emplace(newState->second.value,rules[value].left.id,reducers[value]->combine(reduced_variables));
                break;
            }
            case Action::Accept:
            {
                result = variable;
                return true;
            }
        }
    }
    return false;
}

std::shared_ptr<Variable> LRLogic::calculate(const std::string &s) {
    parse(s);
    return result;
}

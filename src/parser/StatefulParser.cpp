//
// Created by ramizouari on 20/10/23.
//
#include "parser/StatefulParser.h"

using namespace parser;

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

StatefulLR0ParserBuilder::StatefulLR0ParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds)
{

}

StatefulSLRParserBuilder::StatefulSLRParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds)
{

}

StatefulLALRParserBuilder::StatefulLALRParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds)
{

}

StatefulLRParserBuilder::StatefulLRParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds){

}

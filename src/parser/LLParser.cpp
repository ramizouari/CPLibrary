//
// Created by ramizouari on 26/05/22.
//

#include "parser/LLParser.h"

namespace cp::parser
{
    bool LLParser::parse(const std::string &S)
    {

    }

    LLParser &LLParser::build() &{
        buildFirst();
        buildFollow();
        isTerminal.resize(symbols.size(),true);

        for(const auto & rule : rules)
            isTerminal[rule.left.id]=false;
        for(int i=0;i<rules.size();i++) for(auto S : rules[i].right)
        {
            auto& rule=rules[i];
            for (auto terminalId: firstIds[S.id])
            {
                if (terminalId == SpecialCharacter::Epsilon)
                    continue;
                actionMap.emplace(std::make_pair(S.id, terminalId), Action(Action::Expand, i));
            }
            if(!firstIds[S.id].count(SpecialCharacter::Epsilon))
                break;
        }
    }

    LLFamily::Action::Action(LLFamily::Action::Type type, std::uint64_t value) : type(type), value(value) {}

} // parser
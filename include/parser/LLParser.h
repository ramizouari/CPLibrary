//
// Created by ramizouari on 26/05/22.
//

#ifndef UTF8_LLPARSER_H
#define UTF8_LLPARSER_H
#include "Grammar.h"
namespace cp
{
    namespace parser {

        namespace LLFamily
        {
            struct Action
            {
                enum Type
                {
                    Expand,
                    Accept
                };
                Type type;
                std::uint64_t value;
                Action(Type type,std::uint64_t value);
            };
        }

        class LLParser final: protected GrammarWithFirstFollow,public StringParser<char>
        {
            using Action=LLFamily::Action;
            std::vector<bool> isTerminal;
            std::unordered_map<std::pair<std::uint64_t,std::uint64_t>,Action> actionMap;

        public:
            bool parse(const std::string &s) override;
            LLParser& build() &;
        };

    } // parser
}

#endif //UTF8_LLPARSER_H

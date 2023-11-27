//
// Created by ramizouari on 23/05/22.
//

#include <regex>
#include "parser/Grammar.h"

namespace cp
{
    std::uint64_t Symbol::initHash() {
        return hasher(id);
    }

    Symbol::Symbol(std::uint64_t id): id(id)
    {
        initHash();
    }

    Rule::Rule(Symbol left, const std::vector<Symbol> &right):left(left),right(right)
    {
        initHash();
    }

    std::uint64_t Rule::initHash()  {
        std::uint64_t hash=left.hash;
        hash^=hasher(right);
        return hash;
    }

    Rule::Rule(Symbol left, std::vector<Symbol> &&right) : left(left), right(std::move(right))
    {
        initHash();
    }

    void Grammar::addRule(const Rule &rule) {
        if(rule.left.id>=ruleIdsBySymbolId.size())
            ruleIdsBySymbolId.resize(2*rule.left.id+1);
        ruleIdsBySymbolId[rule.left.id].push_back(rules.size());
        rules.push_back(rule);
    }

    const std::vector<Rule> &Grammar::getRules() const {
        return rules;
    }

    void Grammar::addRule(Rule &&rule) {
        if(rule.left.id>=ruleIdsBySymbolId.size())
            ruleIdsBySymbolId.resize(2*rule.left.id+1);
        ruleIdsBySymbolId[rule.left.id].push_back(rules.size());
        rules.push_back(std::move(rule));
    }

    void Grammar::setAxiom(const Symbol & _axiom) {
        this->axiom=_axiom;
    }


    void GrammarWithFirstFollow::buildFirst() {
        firstIds.resize(symbols.size());
        bool insertion;
        /*
         * Suppose every symbol is terminal, so that each symbol S has the singleton {S} as symbol set.
         * */
        for(int i=0;i<firstIds.size();i++)
            firstIds[i]={ static_cast<std::uint64_t>(i)};
        /*
         * For every non-terminal symbol, set the empty set as symbol its first set.
         * */
        for(auto &rule:getRules())
            firstIds[rule.left.id]={};

        do
        {
            insertion=false;
            for(auto & rule : getRules())
            {
                auto &leftSymbol=symbols[rule.left.id];
                auto &leftFirstIds=firstIds[leftSymbol.id];
                int i;
                for(i=0;i<rule.right.size();i++)
                {
                    auto &rightSymbol = rule.right[i];
                    auto &rightFirstIds = this->firstIds[rightSymbol.id];
                    for(auto &firstId : rightFirstIds)
                        if(firstId != SpecialCharacter::Epsilon && leftFirstIds.insert(firstId).second)
                            insertion=true;
                    if(!rightFirstIds.contains(SpecialCharacter::Epsilon))
                        break;
                }
                if(i==rule.right.size() && leftFirstIds.insert(SpecialCharacter::Epsilon).second)
                    insertion=true;
            }
        } while(insertion);
    }

    void GrammarWithFirstFollow::buildFollow()
    {
        followIds.resize(symbols.size());
        bool insertion;
        followIds[axiom.id].insert(SpecialCharacter::EndOfFile);
        do
        {
            insertion=false;
            for(auto & rule : getRules())
            {
                auto &leftSymbol=symbols[rule.left.id];
                auto &leftFirstIds=firstIds[leftSymbol.id];
                for(int i=0;i<rule.right.size();i++)
                {
                    bool allRightsHaveEpsilon = true;
                    auto &S1 = rule.right[i];
                    for (int j = i + 1; j < rule.right.size(); j++)
                    {
                        auto &S2 = rule.right[j];
                        for(auto &id:firstIds[S2.id]) if(id != SpecialCharacter::Epsilon)
                            {
                                auto [_,inserted]=followIds[S1.id].insert(id);
                                insertion=insertion||inserted;
                            }
                        if (!firstIds[S2.id].contains(SpecialCharacter::Epsilon))
                        {
                            allRightsHaveEpsilon = false;
                            break;
                        }
                    }
                    if(allRightsHaveEpsilon)
                    {
                        for(auto &id:followIds[leftSymbol.id])
                        {
                            auto [_,inserted]=followIds[S1.id].insert(id);
                            insertion=insertion||inserted;
                        }
                    }
                }
            }
        } while(insertion);
    }

}
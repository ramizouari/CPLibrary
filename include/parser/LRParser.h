//
// Created by ramizouari on 26/05/22.
//

#ifndef UTF8_LRPARSER_H
#define UTF8_LRPARSER_H
#include "Grammar.h"
/*
 * This header contains the implementations of LR(1) family parsers.
 * The implementation is based on the theory of the LR(1) parser.
 *
 * The supported parsers are:
 * - SLR(1)
 * - LALR(1)
 * - LR(1)
 * */


namespace parser {

    struct LR0Item : public Hashable
    {
        std::uint64_t ruleId;
        std::uint64_t dot;
        explicit LR0Item(std::uint64_t ruleId,std::uint64_t dot);
        bool operator==(const LR0Item & other) const = default;
        auto operator<=>(const LR0Item & other) const = default;
        std::uint64_t initHash();
    private:
        inline static constexpr std::hash<std::uint64_t> hasher{};
    };

    struct LR1Item : public Hashable
    {
        std::uint64_t ruleId;
        std::uint64_t dot;
        std::uint64_t lookahead;
        explicit LR1Item(std::uint64_t ruleId,std::uint64_t dot,std::uint64_t lookahead);
        bool operator==(const LR1Item & other) const = default;
        auto operator<=>(const LR1Item & other) const = default;
        std::uint64_t initHash();
    private:
        inline static constexpr std::hash<std::uint64_t> hasher{};
    };

    namespace LRFamily
    {
        struct Action
        {
            enum Type
            {
                Shift,
                Reduce,
                Accept,
                Error
            };
            Type type;
            std::uint64_t value;
            Action(Type type,std::uint64_t value);
        };
    }

    class SLRParser final: protected  GrammarWithFirstFollow, public StringParser<char>
    {

        using Action=LRFamily::Action;
        std::unordered_map<std::unordered_set<LR0Item>,int> lr0ItemsIds;
        std::vector<std::unordered_set<LR0Item>> lr0Items;
        std::vector<bool> isTerminal;
        std::unordered_set<LR0Item> closure(const std::unordered_set<LR0Item> &item);
        std::unordered_map<std::pair<std::uint64_t,std::uint64_t>,Action> gotoIds;
    public:
        bool parse(const std::string &s) override;
        SLRParser& build() &;
    };

    class LRParser final: protected  GrammarWithFirstFollow, public StringParser<char>
    {

        using Action=LRFamily::Action;
        std::unordered_map<std::unordered_set<LR1Item>,int> lr1ItemsIds;
        std::vector<std::unordered_set<LR1Item>> lr1Items;
        std::vector<bool> isTerminal;
        std::unordered_set<LR1Item> closure(const std::unordered_set<LR1Item> &item);
        std::unordered_map<std::pair<std::uint64_t,std::uint64_t>,Action> gotoIds;
    public:
        bool parse(const std::string &s) override;
        LRParser& build() &;
    };

    class LALRParser final: protected  GrammarWithFirstFollow, public StringParser<char>
    {

        using Action=LRFamily::Action;
        std::unordered_map<std::unordered_set<LR0Item>,int> lr0ItemsIds;
        std::vector<std::unordered_set<LR0Item>> lr0Items;
        std::vector<bool> isTerminal;
        std::unordered_set<LR0Item> closure(const std::unordered_set<LR0Item> &item);
        std::unordered_map<std::pair<std::uint64_t,std::uint64_t>,Action> gotoIds;
    public:
        bool parse(const std::string &s) override;
        SLRParser& build() &;
    };


} // parser

#endif //UTF8_LRPARSER_H

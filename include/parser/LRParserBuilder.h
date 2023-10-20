//
// Created by ramizouari on 26/05/22.
//

#ifndef UTF8_LRPARSER_H
#define UTF8_LRPARSER_H
#include "Grammar.h"
#include "StatefulParser.h"
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

    struct LALR1Item : public Hashable
    {
        std::uint64_t ruleId;
        std::uint64_t dot;
        mutable std::unordered_set<std::uint64_t> lookahead;
        explicit LALR1Item(std::uint64_t ruleId,std::uint64_t dot,std::unordered_set<std::uint64_t>&& lookahead);
        explicit LALR1Item(std::uint64_t ruleId,std::uint64_t dot,const std::unordered_set<std::uint64_t>& lookahead);
        bool operator==(const LALR1Item & other) const;
        std::uint64_t initHash();
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

        struct LALR1ItemContext
        {
            std::uint64_t id;
            std::unordered_set<std::uint64_t> lookaheads;
            explicit LALR1ItemContext(std::uint64_t id,std::unordered_set<std::uint64_t> &&lookaheads);
            explicit LALR1ItemContext(std::uint64_t id,const std::unordered_set<std::uint64_t> &lookaheads);

        };
    }

    class ShiftReduceParser :virtual public StringParser<char>
    {
    public:
        using Action=LRFamily::Action;
        void printTable(std::ostream &H) const;
        bool parse(const std::string &s) override;
        //bool parse(std::uint64_t symbolId);
    protected:
        std::unordered_map<std::pair<std::uint64_t,std::uint64_t>,Action> gotoIds;
    };

    class LRFamilyBuilder
    {
        virtual LRFamilyBuilder& build() & = 0;
    };

    class LR0ParserBuilder: public ShiftReduceParser, public LRFamilyBuilder
    {

        std::unordered_map<std::unordered_set<LR0Item>,int> lr0ItemsIds;
        std::vector<std::unordered_set<LR0Item>> lr0Items;
        std::vector<bool> isTerminal;
        std::unordered_set<LR0Item> closure(const std::unordered_set<LR0Item> &item);
    public:
        LR0ParserBuilder& build() &;
    };

    class SLRParserBuilder: protected  GrammarWithFirstFollow, public ShiftReduceParser, public LRFamilyBuilder
    {

        using Action=LRFamily::Action;
        std::unordered_map<std::unordered_set<LR0Item>,int> lr0ItemsIds;
        std::vector<std::unordered_set<LR0Item>> lr0Items;
        std::vector<bool> isTerminal;
        std::unordered_set<LR0Item> closure(const std::unordered_set<LR0Item> &item);
    public:
        SLRParserBuilder& build() &;
    };

    class LRParserBuilder: virtual protected  GrammarWithFirstFollow, virtual public ShiftReduceParser, public LRFamilyBuilder
    {

        using Action=LRFamily::Action;
        std::unordered_map<std::unordered_set<LR1Item>,int> lr1ItemsIds;
        std::vector<std::unordered_set<LR1Item>> lr1Items;
        std::vector<bool> isTerminal;
        std::unordered_set<LR1Item> closure(const std::unordered_set<LR1Item> &item);
    public:
        //bool parse(const std::string &s) override;
        LRParserBuilder& build() &;
    };

    class LALRParserBuilder: virtual protected  GrammarWithFirstFollow, virtual public ShiftReduceParser, public LRFamilyBuilder
    {

        using Action=LRFamily::Action;
        using LRLALR1ItemContext=LRFamily::LALR1ItemContext;
        std::unordered_map<std::unordered_set<LALR1Item>,std::uint64_t> lr1ItemsIds;
        std::vector<std::unordered_set<LALR1Item>> lr1Items;
        std::vector<bool> isTerminal;
        std::unordered_set<LALR1Item> closure(const std::unordered_set<LALR1Item> &item);
    public:
        LALRParserBuilder& build() &;
    };

    class ArithmeticParser : public LRParserBuilder
    {

    public:
    };

    class StatefulShiftReduceParser : virtual public StatefulStringParser<char>
    {
    public:
        using Action=LRFamily::Action;
        std::shared_ptr<Variable> evaluate(const std::string &s) override;
        using IdsMapType = std::unordered_map<std::pair<std::uint64_t,std::uint64_t>,Action>;
        explicit StatefulShiftReduceParser(IdsMapType &gotoIds);
    protected:
        IdsMapType& gotoIds;
    };

    struct StatefulLR0ParserBuilder: public LR0ParserBuilder , public StatefulShiftReduceParser
    {
        StatefulLR0ParserBuilder();
    };

    struct StatefulSLRParserBuilder: public SLRParserBuilder , public StatefulShiftReduceParser
    {
        StatefulSLRParserBuilder();
    };

    struct StatefulLALRParserBuilder: virtual public LALRParserBuilder , virtual public StatefulShiftReduceParser
    {
        StatefulLALRParserBuilder();
    };

    struct StatefulLRParserBuilder: virtual public LRParserBuilder , virtual public StatefulShiftReduceParser
    {
        StatefulLRParserBuilder();
    };

} // parser

#endif //UTF8_LRPARSER_H

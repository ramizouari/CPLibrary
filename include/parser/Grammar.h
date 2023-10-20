//
// Created by ramizouari on 23/05/22.
//

#ifndef UTF8_GRAMMAR_H
#define UTF8_GRAMMAR_H
#include <unordered_map>
#include <unordered_set>
#include <cstdint>
#include <vector>
#include <span>
#include <string>
#include <numeric>
#include <regex>
#include <iostream>

struct Hashable
{
    std::uint64_t hash=0;
    bool operator==(const Hashable & other) const=default;
    auto operator<=>(const Hashable & other) const= default;
protected:
};

struct Symbol : public Hashable
{
    std::uint64_t id;
    explicit Symbol(std::uint64_t id=0);
    bool operator==(const Symbol & other) const = default;
    auto operator<=>(const Symbol & other) const = default;
    std::uint64_t initHash();
private:
    inline static constexpr std::hash<std::uint64_t> hasher{};
};


template<typename T>
concept isHashable =  std::is_base_of_v<Hashable, std::remove_cvref_t<T>>;

template<isHashable H >
struct std::hash<H>
{
    std::size_t operator()(const H & h) const
    {
        return h.hash;
    }
};

template<typename T>
struct std::hash<std::vector<T>>
{
    inline static constexpr std::hash<T> hasher{};
    std::size_t operator()(const std::vector<T> & v) const
    {
        std::size_t hash=0;
        for(auto & e : v)
            hash^=hasher(e);
        return hash;
    }
};

template<typename A,typename B>
struct std::hash<std::pair<A,B>>
{
    inline static constexpr std::hash<A> hasherA{};
    inline static constexpr std::hash<B> hasherB{};
    std::size_t operator()(const std::pair<A,B> & v) const
    {
        return hasherA(v.first)^hasherB(v.second);
    }
};

template<typename T>
struct std::hash<std::unordered_set<T>>
{
    inline static constexpr std::hash<T> hasher{};
    std::size_t operator()(const std::unordered_set<T> & v) const
    {
        std::size_t hash=0;
        for(auto & e : v)
            hash^=hasher(e);
        return hash;
    }
};

struct Rule : public Hashable
{
    Symbol left;
    std::vector<Symbol> right;
    Rule(Symbol left, const std::vector<Symbol>& right);
    Rule(Symbol left, std::vector<Symbol>&& right);
    bool operator==(const Rule & other) const = default;
    auto operator<=>(const Rule & other) const = default;
    std::uint64_t initHash();
private:
    inline static constexpr std::hash<std::vector<Symbol>> hasher{};

};

class Grammar
{
protected:
    std::vector<Rule>  rules;
    Symbol axiom;
    std::vector<Symbol> symbols;
    std::vector<std::vector<std::uint64_t>> ruleIdsBySymbolId;
public:
    void addRule(const Rule &rule);
    void addRule(Rule &&rule);
    void setAxiom(const Symbol &_axiom);
    virtual ~Grammar() = default;

    template<typename Left, typename... Right>
    void emplaceRule(Left &&left,Right &&... right)
    {
        rules.template emplace_back(std::forward(left),{std::forward(right)...});
    }
    template<typename Left,typename Right>
    void addRule(Left &&left,Right &&right)
    {
        rules.template emplace_back(std::forward(left),std::forward(right));
    }
    [[nodiscard]] const std::vector<Rule> & getRules() const;
};

enum SpecialCharacter : std::uint64_t
{
    Epsilon=std::numeric_limits<std::uint64_t>::max()-2,
    EndOfString,
    Empty=Epsilon,
    EndOfFile=EndOfString
};

class GrammarWithFirstFollow : virtual public Grammar
{
protected:
    std::vector<std::unordered_set<std::uint64_t>> firstIds,followIds;
    virtual void buildFirst();
    virtual void buildFollow();
};



template<typename CharType,typename Traits=std::char_traits<CharType>>
class StringParser : virtual protected Grammar
{


protected:
    using stringType=std::basic_string<CharType,Traits>;
    static std::vector<stringType> split(const stringType &S,const stringType &regex)
    {
        std::basic_regex<CharType,std::regex_traits<CharType>> rgx(regex);
        std::regex_token_iterator<typename stringType::const_iterator> iter(S.begin(),
                                        S.end(),
                                        rgx,
                                        -1),
                                        end;
        std::vector<stringType> result(iter, end);
        if(result.front().empty())
            return {};
        return result;
    }
    std::unordered_map<Symbol,std::uint64_t> symbolId;
    std::unordered_map<stringType,std::uint64_t> symbolMap;
    Symbol& addSymbol(const stringType & s)
    {
        auto it=symbolMap.find(s);
        if(it==symbolMap.end())
        {
            Symbol S{symbolId.size()};
            symbols.push_back(S);
            symbolMap[s]=symbolId.size();
            symbolId[S]=S.id;
            return symbols.back();
        }
        else return symbols[it->second];
    }

public:

    void addRuleList(const stringType &left,const std::vector<stringType> &right)
    {
        Symbol  L=addSymbol(left);
        std::vector<Symbol> R;
        for(const auto & s : right)
            R.push_back(addSymbol(s));
        Rule rule{L,R};
        Grammar::addRule(std::move(rule));
    }

    template<typename Left,typename... Right>
    void addRuleList(Left && left, Right && ... right)
    {
        auto S=addSymbol(std::forward<Left>(left));
        std::vector<Symbol> rightVec{addSymbol(std::forward<Right>(right))...};
        Rule rule(S,rightVec);
        Grammar::addRule(std::move(rule));
    }

    virtual void addRule(const std::basic_string<CharType,Traits> &line)
    {
        std::basic_regex<CharType,std::regex_traits<CharType>> ruleRegex(R"(\s*([^\s]+)\s*(?:->|:|:==)\s*((?:[^\s]\s+)*(?:[^\s]+)|)\s*)");
        std::match_results<typename stringType::const_iterator> ruleMatch;

        if(std::regex_match(line,ruleMatch,ruleRegex))
        {
            std::vector<Symbol>  rightSymbols;
            Symbol leftSymbol=addSymbol(ruleMatch[1]);
            std::cout << ruleMatch[2] << std::endl;
            for(const auto &s : split(ruleMatch[2],R"(\s+)"))
                 rightSymbols.push_back(addSymbol(s));
            Rule rule(leftSymbol,rightSymbols);
            Grammar::addRule(std::move(rule));
        }
        else
            throw std::runtime_error("Invalid rule");
    }

    //virtual void addRule(std::basic_string<CharType,Traits> &&line) = 0;
    virtual bool parse(const std::basic_string<CharType,Traits> &line)=0;
};





#endif //UTF8_GRAMMAR_H

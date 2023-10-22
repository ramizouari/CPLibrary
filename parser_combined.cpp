//
// Created by ramizouari on 21/10/23.
//
#include <utility>
#include <variant>
#include <list>
#include <valarray>
//
// Created by ramizouari on 26/05/22.
//

#ifndef UTF8_LRPARSER_H
#define UTF8_LRPARSER_H
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
//
// Created by ramizouari on 20/10/23.
//

#ifndef CPLIBRARY_STATEFULPARSER_H
#define CPLIBRARY_STATEFULPARSER_H

namespace parser
{

    struct Variable
    {
        virtual ~Variable() = default;
    };

    struct VariableReducer
    {
        virtual ~VariableReducer() = default;
        virtual std::shared_ptr<Variable> combine(const std::vector<std::shared_ptr<Variable>>& symbols) const = 0;
    };

    namespace reducers
    {
        struct NullReducer final: public VariableReducer
        {
            [[nodiscard]] std::shared_ptr<Variable> combine(const std::vector<std::shared_ptr<Variable>>& symbols) const
            {
                return nullptr;
            }
            inline static std::shared_ptr<VariableReducer> instance = std::make_shared<NullReducer>();
        };

        struct Projection : public VariableReducer
        {
            std::int64_t index;

            [[nodiscard]] std::shared_ptr<Variable> combine(const std::vector<std::shared_ptr<Variable>>& symbols) const
            {
                return symbols.at(index);
            }
            explicit Projection(std::int64_t index) : index(index){}
        };

        struct Identity : public Projection
        {
            Identity(): Projection(0){}
            inline static std::shared_ptr<VariableReducer> instance = std::make_shared<Identity>();
        };

        struct IdentityOrNull : public VariableReducer
        {
            [[nodiscard]] std::shared_ptr<Variable> combine(const std::vector<std::shared_ptr<Variable>>& symbols) const
            {
                return symbols.empty() ? nullptr : symbols.front();
            }
            inline static std::shared_ptr<VariableReducer> instance = std::make_shared<IdentityOrNull>();
        };
    }

    template<typename CharType,typename Traits=std::char_traits<CharType>>
    class StatefulStringParser : virtual protected StringParser<CharType,Traits>
    {
    protected:
        std::vector<std::shared_ptr<VariableReducer>> reducers;
    public:
        using stringType=StringParser<CharType,Traits>::stringType;
        void addRuleList(std::shared_ptr<VariableReducer> reducer,const stringType &left,const std::vector<stringType> &right)
        {
            StringParser<CharType,Traits>::addRule(left,right);
            reducers.push_back(reducer);
        }

        void addRuleList(const stringType &left,const std::vector<stringType> &right)
        {
            addRuleList(reducers::Identity::instance,left,right);
        }



        virtual void addRule(std::shared_ptr<VariableReducer> reducer, const std::basic_string<CharType,Traits> &line)
        {
            StringParser<CharType,Traits>::addRule(line);
            reducers.push_back(reducer);
        }

        template<typename Left,typename... Right>
        void addRuleList(std::shared_ptr<VariableReducer> reducer,Left && left, Right && ... right)
        {
            StringParser<CharType,Traits>::addRuleList(std::forward<Left>(left),std::forward<Right>(right)...);
            reducers.push_back(reducer);
        }


        void addRule(const std::basic_string<CharType,Traits> &line) override
        {
            addRule(reducers::IdentityOrNull::instance,line);
        }

        //virtual void addRule(std::basic_string<CharType,Traits> &&line) = 0;
        virtual std::shared_ptr<Variable> evaluate(const std::basic_string<CharType,Traits> &line)=0;
    };
}

#endif //CPLIBRARY_STATEFULPARSER_H
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
#define NO_ASSIGNMENT ""

struct Function : public parser::VariableReducer
{
    size_t parameters;
};


template<typename T>
struct Type : public parser::Variable
{
    T value;
    explicit Type(T value):value(std::move(value)){}
    Type() = default;
};
using StringValue = Type<std::string>;
using IntValue = Type<int>;

struct CharacterGenerator : public parser::VariableReducer
{
    char c;
    explicit CharacterGenerator(char c):c(c){}
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<StringValue>(std::string(1,c));
    }
};

struct NameConcat : public parser::VariableReducer
{
    char suffix;
    explicit NameConcat(char suffix): suffix(suffix){}

    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto name = std::dynamic_pointer_cast<StringValue>(variables.at(0));
        name->value.push_back(suffix);
        return name;
    }
};

struct ProgramState : public parser::Variable
{
    std::map<std::string, int> variables;
};

struct LastProjection : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return variables.back();
    }
};


struct ListNames : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        using Vector = Type<std::vector<std::string>>;
        auto& top = std::dynamic_pointer_cast<Vector>(v.front())->value;
        auto& last= std::dynamic_pointer_cast<StringValue>(v.back())->value;
        top.push_back(last);
        return v.front();
    }
};

template<typename T>
struct DefaultGenerator : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        return std::make_shared<T>();
    }
};


template<typename T>
struct ToList : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        using Vector = Type<std::vector<T>>;
        auto x= std::make_shared<Vector>();
        x->value.push_back(std::dynamic_pointer_cast<Type<T>>(v.front())->value);
        return x;
    }
};

class ExecutionGraph;
class ExpressionTree;

struct FunctionDefinition : public parser::Variable
{
    std::shared_ptr<StringValue> fn_name;
    std::shared_ptr<Type<std::vector<std::string>>> formal_parameters;
    std::shared_ptr<ExecutionGraph> graph;
    int evaluate(class  Context & context);
};

struct FunctionCall : public parser::Variable
{
    std::string name;
    std::vector<std::shared_ptr<ExpressionTree>> parameters;
};

struct FnDefGenerator : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<FunctionDefinition>();
        x->fn_name = std::dynamic_pointer_cast<StringValue>(v.at(0));
        x->formal_parameters = std::dynamic_pointer_cast<Type<std::vector<std::string>>>(v.at(1));
        x->graph = std::dynamic_pointer_cast<ExecutionGraph>(v.at(2));
        return x;
    }
};

class ExpressionTree;
struct Context
{
    std::map<std::string, int> variables;
    std::map<std::string, std::shared_ptr<FunctionDefinition>> functions;
};

class FunctionCall;

struct ExpressionNode
{
    enum Operation
    {
        Add,
        Sub,
        Mul,
        Div,
        Mod,
        Concat,
        Number,
        Variable,
        FunCall
    };
    Operation op;
    std::shared_ptr<ExpressionNode> lhs,rhs;
    std::variant<std::string, int,FunctionCall> value;
    int evaluate(Context &context);
};

struct ExpressionTree : public parser::Variable
{
    std::shared_ptr<ExpressionNode> root;
    int evaluate(Context &context)
    {
        return root->evaluate(context);
    }
};

int ExpressionNode::evaluate(Context &context)
{
    switch(op)
    {
        case Add:
            return lhs->evaluate(context)+rhs->evaluate(context);
        case Sub:
            return lhs->evaluate(context)-rhs->evaluate(context);
        case Mul:
            return lhs->evaluate(context)*rhs->evaluate(context);
        case Div:
            return lhs->evaluate(context)/rhs->evaluate(context);
        case Mod:
            return lhs->evaluate(context)%rhs->evaluate(context);
        case Concat:
            return lhs->evaluate(context)*std::pow(10,rhs->evaluate(context));
        case Number:
            return std::get<int>(value);
        case Variable:
            return context.variables.at(std::get<std::string>(value));
        case FunCall:
        {
            Context new_context;
            new_context.functions = context.functions;
            auto fn_call = std::get<FunctionCall>(value);
            auto fn = new_context.functions.at(fn_call.name);
            for(int i=0;i<fn->formal_parameters->value.size();i++)
                new_context.variables[fn->formal_parameters->value[i]] = fn_call.parameters[i]->evaluate(context);
            return fn->evaluate(new_context);
        }
    }
}

struct CallParams : public parser::Variable
{
    std::vector<std::shared_ptr<ExpressionTree>> parameters;
};

struct ConcatenateCallParams : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::dynamic_pointer_cast<CallParams>(v.at(0));
        x->parameters.push_back(std::dynamic_pointer_cast<ExpressionTree>(v.at(2)));
        return x;
    }
};

struct Statement : public parser::Variable
{
    std::string lhs;
    std::shared_ptr<ExpressionTree> rhs;
    int evaluate(Context &context)
    {
        return context.variables[lhs] = rhs->evaluate(context);
    }
};


struct ExecutionGraph : public parser::Variable
{
    std::vector<std::shared_ptr<Statement>> statements;
};

int FunctionDefinition::evaluate(Context &context)
{
    int result;
    for(auto statement : graph->statements)
        result=statement->evaluate(context);
    return result;
}

using Number = Type<int>;

struct NumberGenerator : public parser::VariableReducer
{
    std::int64_t v;
public:
    NumberGenerator(std::int64_t v):v(v){}
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(v);
    }
};

struct NumberConcatenation : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(10*std::dynamic_pointer_cast<Number>(variables.at(0))->value + std::dynamic_pointer_cast<Number>(variables.at(1))->value);
    }
};

struct ExpressionTreeBuilder : public parser::VariableReducer
{
    ExpressionNode::Operation op;
    explicit ExpressionTreeBuilder(ExpressionNode::Operation op):op(op){}

    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto x= std::make_shared<ExpressionTree>();
        x->root = std::make_shared<ExpressionNode>();
        x->root->op = op;
        x->root->lhs = std::dynamic_pointer_cast<ExpressionTree>(variables.at(0))->root;
        x->root->rhs = std::dynamic_pointer_cast<ExpressionTree>(variables.at(2))->root;
        return x;
    }
};

struct FunctionCallBuilder : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto tree = std::make_shared<ExpressionTree>();
        FunctionCall fn_call;
        fn_call.name = std::dynamic_pointer_cast<StringValue>(v.at(0))->value;
        if(v.size()==4)
            fn_call.parameters = std::dynamic_pointer_cast<CallParams>(v.at(2))->parameters;
        else if(v.size()!=3)
            throw std::runtime_error("invalid function call");
        tree->root = std::make_shared<ExpressionNode>();
        tree->root->lhs=nullptr;
        tree->root->rhs=nullptr;
        tree->root->op = ExpressionNode::FunCall;
        tree->root->value.emplace<FunctionCall>(std::move(fn_call));
        return tree;
    }
};

struct NumberToLeaf : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<ExpressionTree>();
        x->root = std::make_shared<ExpressionNode>();
        x->root->op = ExpressionNode::Number;
        x->root->value = std::dynamic_pointer_cast<Number>(v.at(0))->value;
        return x;
    }
};

struct VariableToLeaf : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<ExpressionTree>();
        x->root = std::make_shared<ExpressionNode>();
        x->root->op = ExpressionNode::Variable;
        x->root->value = std::dynamic_pointer_cast<StringValue>(v.at(0))->value;
        return x;
    }
};

struct StatementBuilder : public parser::VariableReducer
{
    std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<Statement>();
        x->lhs = std::dynamic_pointer_cast<StringValue>(v.at(0))->value;
        x->rhs = std::dynamic_pointer_cast<ExpressionTree>(v.at(2));
        return x;
    }
};

struct StatementConcat : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::dynamic_pointer_cast<ExecutionGraph>(v.at(0));
        x->statements.push_back(std::dynamic_pointer_cast<Statement>(v.at(2)));
        return x;
    }
};

struct FunctionsBlock : public parser::Variable
{
    std::map<std::string, std::shared_ptr<FunctionDefinition>> functions;
};

struct FunctionsBlockConcat : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::dynamic_pointer_cast<FunctionsBlock>(v.at(0));
        auto fn_def=std::dynamic_pointer_cast<FunctionDefinition>(v.at(2));
        x->functions.insert(std::make_pair(fn_def->fn_name->value, fn_def));
        return x;
    }
};

struct FunctionsBlockBuilder : public parser::VariableReducer
{
    std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<FunctionsBlock>();
        x->functions.insert(std::make_pair(std::dynamic_pointer_cast<FunctionDefinition>(v.at(0))->fn_name->value, std::dynamic_pointer_cast<FunctionDefinition>(v.at(0))));
        return x;
    }
};

struct MainProgram : public parser::Variable
{
    std::shared_ptr<Statement> main;
    std::shared_ptr<Context> context;
    [[nodiscard]] int evaluate()
    {
        return main->evaluate(*context);
    }
};

struct PrepareExecution : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<MainProgram>();
        x->main = std::make_shared<Statement>();
        x->main->lhs=NO_ASSIGNMENT;
        x->main->rhs=std::dynamic_pointer_cast<ExpressionTree>(v.at(2));
        x->context = std::make_shared<Context>();
        x->context->functions = std::move(std::dynamic_pointer_cast<FunctionsBlock>(v.at(0))->functions);
        return x;
    }
};

struct PrepareFunctionBody :public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        std::shared_ptr<ExecutionGraph> x;
        if(v.size()>1)
        {
            x=std::dynamic_pointer_cast<ExecutionGraph>(v.front());
            x->statements.push_back(std::make_shared<Statement>());
            x->statements.back()->lhs=NO_ASSIGNMENT;
            x->statements.back()->rhs=std::dynamic_pointer_cast<ExpressionTree>(v.at(2));
        }
        else
        {
            x=std::make_shared<ExecutionGraph>();
            x->statements.push_back(std::make_shared<Statement>());
            x->statements.back()->lhs=NO_ASSIGNMENT;
            x->statements.back()->rhs=std::dynamic_pointer_cast<ExpressionTree>(v.at(0));
        }

        return x;
    }
};

struct ToListStatement : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<ExecutionGraph>();
        x->statements.push_back(std::dynamic_pointer_cast<Statement>(v.at(0)));
        return x;
    }
};

struct ToCallParams : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<CallParams>();
        x->parameters.push_back(std::dynamic_pointer_cast<ExpressionTree>(v.front()));
        return x;
    }
};

class Tokenizer
{
public:
    std::shared_ptr<parser::StatefulShiftReduceParser> builder;

    std::string tokenize(const std::string &source)
    {
        std::string result;
        result.reserve(source.size());
        int start,end;
        for(start=0;start<source.size() && std::isspace(source[start]);start++);
        for(end=source.size(); end > start && std::isspace(source[end-1]);end--);
        for(int i=start;i< end; i++)
        {
            if(!std::isspace(source[i]) || source[i]=='\n')
                result.push_back(source[i]);
        }
        return result;
    }
    auto evaluate(const std::string & source)
    {
        return builder->evaluate(tokenize(source));
    }
};

int main(int argc, char** argv)
{
    using namespace parser;
    using namespace reducers;
    auto G = std::make_shared<StatefulLRParserBuilder>();
    std::string scope;

    G->addRuleList(std::make_shared<LastProjection>(),"Start","Program");
    G->addRuleList(std::make_shared<PrepareExecution>(),"Program","FunDefBlock","EOL", "Expression");
    G->addRuleList(std::make_shared<FunctionsBlockBuilder>(),"FunDefBlock","FunDef");
    G->addRuleList(std::make_shared<FunctionsBlockConcat>(),"FunDefBlock","FunDefBlock","EOL","FunDef");
    G->addRuleList(reducers::Identity::instance,"EOL","EOL","\n");
    G->addRuleList(reducers::Identity::instance,"EOL","\n");
    G->addRuleList(std::make_shared<FnDefGenerator>(),"FunDef","Name","FunParams","FunBody");
    G->addRuleList(std::make_shared<NameConcat>('f'),"Name","Name","f");
    G->addRuleList(std::make_shared<NameConcat>('g'),"Name","Name","g");
    G->addRuleList(std::make_shared<Projection>(1),"FunParams","(","FunParamGroup",")","EOL");
    G->addRuleList(std::make_shared<DefaultGenerator<Type<std::vector<std::string>>>>(),"FunParams","(",")","EOL");
    G->addRuleList(std::make_shared<ListNames>(),"FunParamGroup","FunParamGroup","Comma","Name");
    G->addRuleList(std::make_shared<ToList<std::string>>(),"FunParamGroup","Name");
    G->addRuleList(std::make_shared<PrepareFunctionBody>(),"FunBody","StatementsBlock","EOL","Expression");
    G->addRuleList(std::make_shared<PrepareFunctionBody>(),"FunBody","Expression");
    G->addRuleList(std::make_shared<ToListStatement>(), "StatementsBlock","Statement");
    G->addRuleList(std::make_shared<StatementConcat>(), "StatementsBlock","StatementsBlock","EOL","Statement");
    G->addRuleList(parser::reducers::Identity::instance, "Comma",",");
    G->addRuleList(std::make_shared<StatementBuilder>(), "Statement","Name","=","Expression");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Add), "Expression","Expression","+","Term");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Sub), "Expression","Expression","-","Term");
    G->addRuleList(parser::reducers::Identity::instance, "Expression","Term");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Mul), "Term","Term","*","Factor");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Div), "Term","Term","/","Factor");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Mod), "Term","Term","%","Factor");
    G->addRuleList(parser::reducers::Identity::instance, "Term","Factor");
    G->addRuleList(std::make_shared<NumberToLeaf>(), "Factor","Number");
    G->addRuleList(std::make_shared<VariableToLeaf>(), "Factor","Name");
    G->addRuleList(std::make_shared<FunctionCallBuilder>(), "Factor","Name","(","CallParams",")");
    G->addRuleList(std::make_shared<FunctionCallBuilder>(), "Factor","Name","(",")");
    G->addRuleList(std::make_shared<ConcatenateCallParams>(), "CallParams","CallParams","Comma","Expression");
    G->addRuleList(std::make_shared<ToCallParams>(), "CallParams","Expression");
    G->addRuleList(std::make_shared<parser::reducers::Projection>(1), "Factor","(","Expression",")");
    G->addRuleList(std::make_shared<NumberConcatenation>(), "Number","Number","Digit");
    G->addRuleList(parser::reducers::Identity::instance, "Number","Digit");
    for(int i=0;i<10;i++)
        G->addRuleList(std::make_shared<NumberGenerator>(i), "Digit",std::to_string(i));
    for(char x='a';x<='z';x++)
        G->addRuleList(std::make_shared<CharacterGenerator>(x),"Name",std::string(1,x));
    for(char x='A';x<='Z';x++)
        G->addRuleList(std::make_shared<CharacterGenerator>(x),"Name",std::string(1,x));
    G->build();
//    G->printTable(std::cout);
    Tokenizer T;
    T.builder = G;
    std::string inputs;
    std::string source_code;
    int n;
    std::cin >> n;
    std::cin.ignore();
    while(n--)
    {
        std::getline(std::cin,inputs);
        source_code.push_back('\n');
        source_code+=inputs;
    }
    auto result = T.evaluate(source_code);
    if(result)
    {
        auto main = std::dynamic_pointer_cast<MainProgram>(result);
        std::cout << main->evaluate() << '\n';
    }
    else
    {
        std::cout << "error\n";
    }
}



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

#include <queue>

namespace parser
{
    std::unordered_set<LALR1Item> LALRParserBuilder::closure(const std::unordered_set<LALR1Item> &items)
    {
        std::unordered_set<LALR1Item> newItems;
        struct context
        {
            std::uint64_t symbolId;
            std::unordered_set<std::uint64_t> firstIds;
        };
        std::stack<context> stack;
        for(const auto& item:items)
        {
            auto [it,inserted] = newItems.insert(item);
            if(!inserted)
                it->lookahead.insert(item.lookahead.begin(),item.lookahead.end());

            if(item.dot < rules[item.ruleId].right.size())
            {
                if(isTerminal[rules[item.ruleId].right[item.dot].id])
                    continue;
                context ctx;
                ctx.symbolId=rules[item.ruleId].right[item.dot].id;
                if(item.dot+1==rules[item.ruleId].right.size())
                    ctx.firstIds.insert(item.lookahead.begin(),item.lookahead.end());
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
                    auto itemIterator=newItems.find(LALR1Item(ruleId, 0, {}));
                    if(itemIterator==newItems.end())
                    {
                        newItems.insert(LALR1Item(ruleId, 0, lookaheads));
                        anyInsertion=true;
                    }
                    else for(auto lookahead:lookaheads) if(!itemIterator->lookahead.contains(lookahead))
                            {
                                itemIterator->lookahead.insert(lookahead);
                                anyInsertion=true;
                            }
                    if (anyInsertion && !rules[ruleId].right.empty() && !isTerminal[rules[ruleId].right.front().id])
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
                            ctx.firstIds.insert(lookaheads.begin(),lookaheads.end());
                        ctx.firstIds.erase(SpecialCharacter::Epsilon);
                    }
                    if(anyInsertion)
                        stack.push(ctx);
                }
            }
        }
        return newItems;
    }

    LALRParserBuilder &LALRParserBuilder::build() &{
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
        LALR1Item I0{rules.size()-1,0,{SpecialCharacter::EndOfString}};
        lr1Items.push_back(std::move(closure({I0})));
        lr1ItemsIds.emplace(lr1Items.back(), lr1Items.size() - 1);
        /*
         * Queue of Item-sets that need to be updated.
         * */
        std::queue<int> updates;
        updates.emplace(0);
        int n=1;
        do
        {
            auto k = updates.front();
            updates.pop();
            std::unordered_map<std::uint64_t,std::unordered_set<LALR1Item>> nextItemsSets;
            auto update=lr1ItemsIds.find(lr1Items[k]);
            lr1Items[k] = update->first;
            for(const auto & item : lr1Items[k])
            {
                if(item.dot<rules[item.ruleId].right.size())
                {
                    auto [it,inserted] = nextItemsSets[rules[item.ruleId].right[item.dot].id].emplace(item.ruleId, item.dot + 1, item.lookahead);
                    if(!inserted)
                        throw std::runtime_error("Debugging: duplicate item in LALR1 closure");
                }
            }
            for(auto &[symbolId,itemsSet]: nextItemsSets)
            {
                itemsSet = std::move(closure(itemsSet));
                auto it=lr1ItemsIds.find(itemsSet);
                if(it == lr1ItemsIds.end())
                {
                    lr1Items.push_back(std::move(itemsSet));
                    std::tie(it,std::ignore)=lr1ItemsIds.emplace(lr1Items.back(), lr1Items.size() - 1);
                    updates.emplace(it->second);
                    n++;
                }
                    //Merge lookaheads
                else
                {
                    bool insertion=false;
                    for (auto &item: it->first)
                    {
                        auto item2 = itemsSet.find(item);
                        if (item2 != itemsSet.end()) for (auto lookahead : item2->lookahead) if (!item.lookahead.contains(lookahead))
                                {
                                    item.lookahead.insert(lookahead);
                                    insertion=true;
                                }
                    }
                    if(insertion)
                        updates.emplace(it->second);
                }
                gotoIds.emplace(std::make_pair(k,symbolId),Action(Action::Shift,it->second));
            }
        } while(!updates.empty());
        for(int i=0;i<n;i++) for(const auto &item: lr1Items[i]) if(rules[item.ruleId].right.size() == item.dot)
                {
                    if (item.ruleId == augmentedRuleId)
                        gotoIds.emplace(std::make_pair(i, SpecialCharacter::EndOfString), Action(Action::Accept, 0));
                    else for(auto lookahead : item.lookahead)
                            gotoIds.emplace(std::make_pair(i, lookahead), Action(Action::Reduce, item.ruleId));
                }
        return *this;
    }


    StatefulLALRParserBuilder::StatefulLALRParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds)
    {

    }
}


namespace parser
{
    std::unordered_set<LR0Item> LR0ParserBuilder::closure(const std::unordered_set<LR0Item> &items) {
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

    LR0ParserBuilder &LR0ParserBuilder::build() &
    {
        Symbol augmentedSymbol=addSymbol("");
        Grammar::addRule(Rule(augmentedSymbol,{axiom}));
        setAxiom(augmentedSymbol);
        auto augmentedRuleId=rules.size()-1;
        ruleIdsBySymbolId.resize(symbols.size());
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
                {
                    for (auto &F: symbols) if (isTerminal[F.id])
                        {
                            if (item.ruleId == augmentedRuleId)
                                gotoIds.emplace(std::make_pair(i, SpecialCharacter::EndOfString),
                                                Action(Action::Accept, 0));
                            else gotoIds.emplace(std::make_pair(i, F.id), Action(Action::Reduce, item.ruleId));
                        }
                    gotoIds.emplace(std::make_pair(i, SpecialCharacter::EndOfString), Action(Action::Reduce, item.ruleId));
                }
        return *this;
    }

    StatefulLR0ParserBuilder::StatefulLR0ParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds)
    {

    }

}


namespace parser {

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


namespace parser
{
    SLRParserBuilder& SLRParserBuilder::build() &
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

    std::unordered_set<LR0Item> SLRParserBuilder::closure(const std::unordered_set<LR0Item> &items)
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

    StatefulSLRParserBuilder::StatefulSLRParserBuilder() : StatefulShiftReduceParser(ShiftReduceParser::gotoIds)
    {

    }


}
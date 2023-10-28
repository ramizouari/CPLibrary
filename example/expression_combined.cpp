#include <utility>
#include <variant>
#include <list>
#include <valarray>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <span>
#include <string>
#include <numeric>
#include <regex>
#include <iostream>
#include <set>

template<typename T>
using VE=std::vector<T>;
template<typename T>
using US=std::unordered_set<T>;
template<typename T,typename V>
using UM=std::unordered_map<T,V>;



struct Hashable
{
    std::uint64_t hash=0;
    bool operator==(const Hashable & other) const=default;
    auto operator<=>(const Hashable & other) const= default;
protected:
};

struct SB : public Hashable
{
    std::uint64_t id;
    explicit SB(std::uint64_t id=0);
    bool operator==(const SB & other) const = default;
    auto operator<=>(const SB & other) const = default;
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
struct std::hash<VE<T>>
{
    inline static constexpr std::hash<T> hasher{};
    std::size_t operator()(const VE<T> & v) const
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
struct std::hash<US<T>>
{
    inline static constexpr std::hash<T> hasher{};
    std::size_t operator()(const US<T> & v) const
    {
        std::size_t hash=0;
        for(auto & e : v)
            hash^=hasher(e);
        return hash;
    }
};

struct Rule : public Hashable
{
    SB left;
    VE<SB> right;
    Rule(SB left, const VE<SB>& right);
    Rule(SB left, VE<SB>&& right);
    bool operator==(const Rule & other) const = default;
    auto operator<=>(const Rule & other) const = default;
    std::uint64_t initHash();
private:
    inline static constexpr std::hash<VE<SB>> hasher{};

};

class Grammar
{
protected:
    VE<Rule>  rules;
    SB axiom;
    VE<SB> SBs;
    VE<VE<std::uint64_t>> ruleIdsBySBId;
public:
    void addRule(const Rule &rule);
    void addRule(Rule &&rule);
    void setAxiom(const SB &_axiom);
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
    [[nodiscard]] const VE<Rule> & getRules() const;
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
    VE<US<std::uint64_t>> firstIds,followIds;
    virtual void buildFirst();
    virtual void buildFollow();
};

template<typename CharType,typename Traits=std::char_traits<CharType>>
class SP : virtual protected Grammar
{


protected:
    using stringType=std::basic_string<CharType,Traits>;
    static VE<stringType> split(const stringType &S,const stringType &regex)
    {
        std::basic_regex<CharType,std::regex_traits<CharType>> rgx(regex);
        std::regex_token_iterator<typename stringType::const_iterator> iter(S.begin(),
                                        S.end(),
                                        rgx,
                                        -1),
                                        end;
        VE<stringType> result(iter, end);
        if(result.front().empty())
            return {};
        return result;
    }
    UM<SB,std::uint64_t> SBId;
    UM<stringType,std::uint64_t> SBMap;
    SB& addSB(const stringType & s)
    {
        auto it=SBMap.find(s);
        if(it==SBMap.end())
        {
            SB S{SBId.size()};
            SBs.push_back(S);
            SBMap[s]=SBId.size();
            SBId[S]=S.id;
            return SBs.back();
        }
        else return SBs[it->second];
    }

public:

    void RL(const stringType &left,const VE<stringType> &right)
    {
        SB  L=addSB(left);
        VE<SB> R;
        for(const auto & s : right)
            R.push_back(addSB(s));
        Rule rule{L,R};
        Grammar::addRule(std::move(rule));
    }

    template<typename Left,typename... Right>
    void RL(Left && left, Right && ... right)
    {
        auto S=addSB(std::forward<Left>(left));
        VE<SB> rightVec{addSB(std::forward<Right>(right))...};
        Rule rule(S,rightVec);
        Grammar::addRule(std::move(rule));
    }

    virtual void addRule(const std::basic_string<CharType,Traits> &line)
    {
        std::basic_regex<CharType,std::regex_traits<CharType>> ruleRegex(R"(\s*([^\s]+)\s*(?:->|:|:==)\s*((?:[^\s]\s+)*(?:[^\s]+)|)\s*)");
        std::match_results<typename stringType::const_iterator> ruleMatch;

        if(std::regex_match(line,ruleMatch,ruleRegex))
        {
            VE<SB>  rightSBs;
            SB leftSB=addSB(ruleMatch[1]);
            std::cout << ruleMatch[2] << std::endl;
            for(const auto &s : split(ruleMatch[2],R"(\s+)"))
                 rightSBs.push_back(addSB(s));
            Rule rule(leftSB,rightSBs);
            Grammar::addRule(std::move(rule));
        }
        else
            throw std::runtime_error("Invalid rule");
    }
    virtual bool parse(const std::basic_string<CharType,Traits> &line)=0;
};

namespace parser
{

    struct Variable
    {
        virtual ~Variable() = default;
    };

    struct VariableReducer
    {
        virtual ~VariableReducer() = default;
        virtual std::shared_ptr<Variable> combine(const VE<std::shared_ptr<Variable>>& SBs) const = 0;
    };

    namespace reducers
    {
        struct NullReducer final: public VariableReducer
        {
            [[nodiscard]] std::shared_ptr<Variable> combine(const VE<std::shared_ptr<Variable>>& SBs) const
            {
                return nullptr;
            }
            inline static std::shared_ptr<VariableReducer> instance = std::make_shared<NullReducer>();
        };

        struct Projection : public VariableReducer
        {
            std::int64_t index;

            [[nodiscard]] std::shared_ptr<Variable> combine(const VE<std::shared_ptr<Variable>>& SBs) const
            {
                return SBs.at(index);
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
            [[nodiscard]] std::shared_ptr<Variable> combine(const VE<std::shared_ptr<Variable>>& SBs) const
            {
                return SBs.empty() ? nullptr : SBs.front();
            }
            inline static std::shared_ptr<VariableReducer> instance = std::make_shared<IdentityOrNull>();
        };
    }

    template<typename CharType,typename Traits=std::char_traits<CharType>>
    class SSP : virtual protected SP<CharType,Traits>
    {
    protected:
        VE<std::shared_ptr<VariableReducer>> reducers;
    public:
        using stringType=SP<CharType,Traits>::stringType;
        void RL(std::shared_ptr<VariableReducer> reducer,const stringType &left,const VE<stringType> &right)
        {
            SP<CharType,Traits>::addRule(left,right);
            reducers.push_back(reducer);
        }

        void RL(const stringType &left,const VE<stringType> &right)
        {
            RL(reducers::Identity::instance,left,right);
        }

        virtual void addRule(std::shared_ptr<VariableReducer> reducer, const std::basic_string<CharType,Traits> &line)
        {
            SP<CharType,Traits>::addRule(line);
            reducers.push_back(reducer);
        }

        template<typename Left,typename... Right>
        void RL(std::shared_ptr<VariableReducer> reducer,Left && left, Right && ... right)
        {
            SP<CharType,Traits>::RL(std::forward<Left>(left),std::forward<Right>(right)...);
            reducers.push_back(reducer);
        }


        void addRule(const std::basic_string<CharType,Traits> &line) override
        {
            addRule(reducers::IdentityOrNull::instance,line);
        }

        virtual std::shared_ptr<Variable> evaluate(const std::basic_string<CharType,Traits> &line)=0;
    };
}

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

    class ShiftReduceParser :virtual public SP<char>
    {
    public:
        using Action=LRFamily::Action;
        void printTable(std::ostream &H) const;
        bool parse(const std::string &s) override;
        //bool parse(std::uint64_t SBId);
    protected:
        UM<std::pair<std::uint64_t,std::uint64_t>,Action> gotoIds;
    };

    class LRFamilyBuilder
    {
        virtual LRFamilyBuilder& build() & = 0;
    };


    class LRParserBuilder: virtual protected  GrammarWithFirstFollow, virtual public ShiftReduceParser, public LRFamilyBuilder
    {

        using Action=LRFamily::Action;
        UM<US<LR1Item>,int> lr1ItemsIds;
        VE<US<LR1Item>> lr1Items;
        VE<bool> isTinal;
        US<LR1Item> closure(const US<LR1Item> &item);
    public:
        //bool parse(const std::string &s) override;
        LRParserBuilder& build() &;
    };

    class SSRP : virtual public SSP<char>
    {
    public:
        using Action=LRFamily::Action;
        std::shared_ptr<Variable> evaluate(const std::string &s) override;
        using IdsMapType = UM<std::pair<std::uint64_t,std::uint64_t>,Action>;
        explicit SSRP(IdsMapType &gotoIds);
    protected:
        IdsMapType& gotoIds;
    };

    struct SLRPB: virtual public LRParserBuilder , virtual public SSRP
    {
        SLRPB();
    };

} // parser

template<typename T>
struct Type : public parser::Variable
{
    T value;
    explicit Type(T value):value(std::move(value)){}
    Type() = default;
};
using StringValue = Type<std::string>;
using IntValue = Type<std::int64_t>;


struct LastProjection : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return variables.back();
    }
};

template<typename T>
struct DefaultGenerator : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> &v) const override
    {
        return std::make_shared<T>();
    }
};


template<typename T>
struct ToList : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> &v) const override
    {
        using Vector = Type<VE<T>>;
        auto x= std::make_shared<Vector>();
        x->value.push_back(std::dynamic_pointer_cast<Type<T>>(v.front())->value);
        return x;
    }
};

struct ENode
{
    enum Operation
    {
        Add,
        Sub,
        Mul,
        Div,
        Mod,
        Pow,
        Concat,
        Number,
    };
    Operation op;
    std::shared_ptr<ENode> lhs,rhs;
    std::int64_t value;
    std::int64_t evaluate();
};

struct ETree : public parser::Variable
{
    std::shared_ptr<ENode> root;
    std::int64_t evaluate()
    {
        return root->evaluate();
    }
};

std::int64_t power(std::int64_t a, std::int64_t n)
{
    if(n==0)
        return 1;
    else if(n==1)
        return a;
    auto s=power(a,n/2);
    return n%2?s*s*a:s*s;
}

std::int64_t ENode::evaluate()
{
    switch(op)
    {
        case Add:
            return lhs->evaluate()+rhs->evaluate();
        case Sub:
            return lhs->evaluate()-rhs->evaluate();
        case Mul:
            return lhs->evaluate()*rhs->evaluate();
        case Div:
            return lhs->evaluate()/rhs->evaluate();
        case Mod:
            return lhs->evaluate()%rhs->evaluate();
        case Concat:
            return lhs->evaluate()*std::pow(10,rhs->evaluate());
        case Pow:
            return power(lhs->evaluate(),rhs->evaluate());
        case Number:
            return value;

    }
}

struct CallParams : public parser::Variable
{
    VE<std::shared_ptr<ETree>> parameters;
};

struct ConcatenateCallParams : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::dynamic_pointer_cast<CallParams>(v.at(0));
        x->parameters.push_back(std::dynamic_pointer_cast<ETree>(v.at(2)));
        return x;
    }
};

using Number = Type<std::int64_t>;

struct NumberGenerator : public parser::VariableReducer
{
    std::int64_t v;
public:
    NumberGenerator(std::int64_t v):v(v){}
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<StringValue>(std::to_string(v));
    }
};

struct NumberConcatenation : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto Z=std::dynamic_pointer_cast<StringValue>(variables.at(0));
        auto &x=Z->value;
        x+=std::dynamic_pointer_cast<StringValue>(variables.at(1))->value;
        return Z;
    }
};

struct LeadingNumberConcatenation : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto &x=std::dynamic_pointer_cast<StringValue>(variables.at(0))->value;
        if(variables.size()==1)
            return std::make_shared<Number>(std::stoll(x));
        auto &y=std::dynamic_pointer_cast<StringValue>(variables.at(1))->value;
        auto z=std::stoll(x+y);
        return std::make_shared<Number>(z);
    }
};

struct ETreeBuilder : public parser::VariableReducer
{
    ENode::Operation op;
    explicit ETreeBuilder(ENode::Operation op):op(op){}

    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto x= std::make_shared<ETree>();
        x->root = std::make_shared<ENode>();
        x->root->op = op;
        x->root->lhs = std::dynamic_pointer_cast<ETree>(variables.at(0))->root;
        x->root->rhs = std::dynamic_pointer_cast<ETree>(variables.at(2))->root;
        return x;
    }
};

struct NumberToLeaf : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const VE<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<ETree>();
        x->root = std::make_shared<ENode>();
        x->root->op = ENode::Number;
        x->root->value = std::dynamic_pointer_cast<Number>(v.at(0))->value;
        return x;
    }
};

std::string tokenize(const std::string &in)
{
    std::string out;
    int balance=0;
    std::set<char> operators={'+','*'};
    for(auto c:in)
    {
        out.push_back(c);
        if(c=='(')
            balance++;
        else if(c==')')
            balance--;
        if(balance==0 && ! operators.count(c))
            out.push_back('#');
    }
}

int main(int argc, char** argv)
{
    using namespace parser;
    using namespace reducers;
    auto G = std::make_shared<SLRPB>();
    std::string scope;
    int K,N;
    std::cin >> K >> N;
    G->RL(std::make_shared<LastProjection>(),"Start","E");
    G->RL(std::make_shared<ETreeBuilder>(ENode::Add), "E","E","+","T");
    G->RL(parser::reducers::Identity::instance, "E","T");
    G->RL(std::make_shared<ETreeBuilder>(ENode::Mul), "T","T","*","Factor");
    G->RL(parser::reducers::Identity::instance, "T","Factor");
    G->RL(std::make_shared<NumberToLeaf>(), "Factor","Number");
    G->RL(std::make_shared<parser::reducers::Projection>(1), "Factor","(","E",")");
    G->RL(std::make_shared<LeadingNumberConcatenation>(), "Number","NonZero","NumberWithZero"); // Number -> Number Digit
    G->RL(std::make_shared<LeadingNumberConcatenation>(), "Number","NonZero"); // Number -> Number Digit
    G->RL(std::make_shared<LeadingNumberConcatenation>(), "Number","Zero"); // Number -> Number Digit
    G->RL(parser::reducers::Identity::instance, "NumberWithZero","Digit");
    G->RL(std::make_shared<NumberConcatenation>(), "NumberWithZero","NumberWithZero","Digit");

    for(int i=0;i<10;i++)
        G->RL(std::make_shared<NumberGenerator>(i), "Digit",std::to_string(i)); // Digit -> i
    for(int i=1;i<10;i++)
        G->RL(std::make_shared<NumberGenerator>(i), "NonZero",std::to_string(i)); // Digit -> i
    G->RL(std::make_shared<NumberGenerator>(0), "Zero","0"); // Digit -> i
    G->build();
    //G->printTable(std::cout);
    std::string E;
    std::cin >> E;
    int counter=0;
    for(int i=0;i<N;i++) for(int j=i+1;j<=N;j++)
    {
        auto sub=E.substr(i,j-i);
        auto result = G->evaluate(sub);
        if(result)
        {
            auto main = std::dynamic_pointer_cast<ETree>(result);
            auto value = main->evaluate();
            counter += value%K==0;
        }
    }
    std::cout << counter << '\n';

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
        SB augmentedSB=addSB("");
        Grammar::addRule(Rule(augmentedSB,{axiom}));
        setAxiom(augmentedSB);
        auto augmentedRuleId=rules.size()-1;
        ruleIdsBySBId.resize(SBs.size());
        buildFirst();
        buildFollow();
        isTinal.resize(SBs.size(),true);

        for(const auto & rule : rules)
            isTinal[rule.left.id]=false;
        LR1Item I0{rules.size()-1,0,{SpecialCharacter::EndOfString}};
        lr1Items.push_back(std::move(closure({I0})));
        lr1ItemsIds.emplace(lr1Items.back(), lr1Items.size() - 1);
        int k=0,n=1;
        do
        {
            UM<std::uint64_t,US<LR1Item>> nextItemsSets;
            for(const auto & item : lr1Items[k])
            {
                if(item.dot<rules[item.ruleId].right.size())
                    nextItemsSets[rules[item.ruleId].right[item.dot].id].emplace(item.ruleId,item.dot+1,item.lookahead);
            }
            for(auto &[SBId,itemsSet]: nextItemsSets)
            {
                itemsSet = std::move(closure(itemsSet));
                auto it=lr1ItemsIds.find(itemsSet);
                if(it == lr1ItemsIds.end())
                {
                    lr1Items.push_back(std::move(itemsSet));
                    std::tie(it,std::ignore)=lr1ItemsIds.emplace(lr1Items.back(), lr1Items.size() - 1);
                    n++;
                }
                gotoIds.emplace(std::make_pair(k,SBId),Action(Action::Shift,it->second));
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

    US<LR1Item> LRParserBuilder::closure(const US<LR1Item> &items) {
        US<LR1Item> newItems;
        struct context
        {
            std::uint64_t SBId;
            US<std::uint64_t> firstIds;
        };
        std::stack<context> stack;
        for(const auto& item:items)
        {
            newItems.insert(item);
            if(item.dot < rules[item.ruleId].right.size())
            {
                if(isTinal[rules[item.ruleId].right[item.dot].id])
                    continue;
                context ctx;
                ctx.SBId=rules[item.ruleId].right[item.dot].id;
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
                for(const auto & ruleId : ruleIdsBySBId[id])
                {
                    context ctx;
                    if(!rules[ruleId].right.empty())
                        ctx.SBId = rules[ruleId].right.front().id;
                    for (auto lookahead: lookaheads)
                    {
                        auto [newItemIterator, inserted] = newItems.insert(LR1Item{ruleId, 0, lookahead});
                        anyInsertion = anyInsertion || inserted;
                        if (inserted && !rules[ruleId].right.empty() && !isTinal[rules[ruleId].right.front().id])
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
            auto SBId=SBMap.at(character);
            auto it=gotoIds.find(std::make_pair(state,SBId));
            if(it==gotoIds.end())
                return false;
            auto [type,value]=it->second;
            switch(type)
            {
                case Action::Shift:
                    i++;
                    stack.emplace(value,SBId);
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

    std::shared_ptr<Variable> SSRP::evaluate(const std::string &S)
    {
        std::stack<std::tuple<std::uint64_t,std::uint64_t,std::shared_ptr<Variable>>> stack;
        stack.emplace(0,-1, nullptr);
        char character[2]{0,0};
        for(int i=0;i<S.size();)
        {
            auto &s=S[i];
            auto [state,ruleId,variable]=stack.top();
            character[0]=s;
            auto SBId=SBMap.at(character);
            auto it=gotoIds.find(std::make_pair(state,SBId));
            if(it==gotoIds.end())
                return nullptr;
            auto [type,value]=it->second;
            switch(type)
            {
                case Action::Shift:
                    i++;
                    stack.emplace(value,SBId,variable);
                    break;
                case Action::Reduce:
                {
                    VE<std::shared_ptr<Variable>> reduced_variables(rules[value].right.size());
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
                    VE<std::shared_ptr<Variable>> reduced_variables(rules[value].right.size());
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

    SSRP::SSRP(SSRP::IdsMapType &gotoIds) : gotoIds(gotoIds) {

    }

    SLRPB::SLRPB() : SSRP(ShiftReduceParser::gotoIds){

    }

} 
std::uint64_t SB::initHash() {
    return hasher(id);
}

SB::SB(std::uint64_t id): id(id)
{
    initHash();
}

Rule::Rule(SB left, const VE<SB> &right):left(left),right(right)
{
    initHash();
}

std::uint64_t Rule::initHash()  {
    std::uint64_t hash=left.hash;
    hash^=hasher(right);
    return hash;
}

Rule::Rule(SB left, VE<SB> &&right) : left(left), right(std::move(right))
{
    initHash();
}

void Grammar::addRule(const Rule &rule) {
    if(rule.left.id>=ruleIdsBySBId.size())
        ruleIdsBySBId.resize(2*rule.left.id+1);
    ruleIdsBySBId[rule.left.id].push_back(rules.size());
    rules.push_back(rule);
}

const VE<Rule> &Grammar::getRules() const {
    return rules;
}

void Grammar::addRule(Rule &&rule) {
    if(rule.left.id>=ruleIdsBySBId.size())
        ruleIdsBySBId.resize(2*rule.left.id+1);
    ruleIdsBySBId[rule.left.id].push_back(rules.size());
    rules.push_back(std::move(rule));
}

void Grammar::setAxiom(const SB & _axiom) {
    this->axiom=_axiom;
}


void GrammarWithFirstFollow::buildFirst() {
    firstIds.resize(SBs.size());
    bool insertion;
    for(int i=0;i<firstIds.size();i++)
        firstIds[i]={ static_cast<std::uint64_t>(i)};
    for(auto &rule:getRules())
        firstIds[rule.left.id]={};
    do
    {
        insertion=false;
        for(auto & rule : getRules())
        {
            auto &leftSB=SBs[rule.left.id];
            auto &leftFirstIds=firstIds[leftSB.id];
            int i;
            for(i=0;i<rule.right.size();i++)
            {
                auto &rightSB = rule.right[i];
                auto &rightFirstIds = this->firstIds[rightSB.id];
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
    followIds.resize(SBs.size());
    bool insertion;
    followIds[axiom.id].insert(SpecialCharacter::EndOfFile);
    do
    {
        insertion=false;
        for(auto & rule : getRules())
        {
            auto &leftSB=SBs[rule.left.id];
            auto &leftFirstIds=firstIds[leftSB.id];
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
                    for(auto &id:followIds[leftSB.id])
                    {
                        auto [_,inserted]=followIds[S1.id].insert(id);
                        insertion=insertion||inserted;
                    }
                }
            }
        }
    } while(insertion);
}

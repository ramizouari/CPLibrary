//
// Created by ramizouari on 21/10/23.
//
#include <utility>
#include <variant>
#include <list>
#include <valarray>
#include "parser/LRParserBuilder.h"
#define NO_ASSIGNMENT ""
#define OUTPUT_CODE "output"
#define OUTPUT_SYMBOL "#"
#define RETURN_CODE "return"
#define RETURN_SYMBOL "@"
#include <set>

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
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return variables.back();
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

struct ExpressionNode
{
    enum Operation
    {
        Add,
        Mul,
        Number,
        Parenthesis
    };
    Operation op;
    std::shared_ptr<ExpressionNode> lhs,rhs;
    bool tainted=false;
    int *counter=nullptr;
    int K;
    std::int64_t value;
    std::int64_t evaluate();
    std::int64_t evaluate(std::int64_t left, Operation op);
    void taint();
};

struct ExpressionTree : public parser::Variable
{
    std::shared_ptr<ExpressionNode> root;
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

std::int64_t ExpressionNode::evaluate()
{
    if(tainted && lhs)
    {
        lhs->counter=counter;
        lhs->K=K;
    }
    std::int64_t result;
    switch(op)
    {
        case Add:
            result= lhs->evaluate()+rhs->evaluate();
        case Mul:
            result= lhs->evaluate()*rhs->evaluate();
        case Number:
            result= value;
        case Parenthesis:
            result= lhs->evaluate();

    }
    return result;
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


using Number = Type<std::int64_t>;

struct NumberGenerator : public parser::VariableReducer
{
    std::int64_t v;
public:
    NumberGenerator(std::int64_t v):v(v){}
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<StringValue>(std::to_string(v));
    }
};

struct NumberConcatenation : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto Z=std::dynamic_pointer_cast<StringValue>(variables.at(0));
        auto &x=Z->value;
        x+=std::dynamic_pointer_cast<StringValue>(variables.at(1))->value;
        return Z;
    }
};

struct LeadingNumberConcatenation : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto &x=std::dynamic_pointer_cast<StringValue>(variables.at(0))->value;
        if(variables.size()==1)
            return std::make_shared<Number>(std::stoll(x));
        auto &y=std::dynamic_pointer_cast<StringValue>(variables.at(1))->value;
        auto z=std::stoll(x+y);
        return std::make_shared<Number>(z);
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

struct ParenthesisBuilder : public parser::VariableReducer
{

    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto x= std::make_shared<ExpressionTree>();
        x->root = std::make_shared<ExpressionNode>();
        x->root->op = ExpressionNode::Parenthesis;
        x->root->lhs = std::dynamic_pointer_cast<ExpressionTree>(variables.at(1))->root;
        return x;
    }
};

struct ExpressionTreeBuilderWithCounter : public parser::VariableReducer
{
    ExpressionNode::Operation op;
    int &counter;
    int K;
    explicit ExpressionTreeBuilderWithCounter(ExpressionNode::Operation op,int &counter,int K):op(op),counter(counter),K(K){}

    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        auto x= std::make_shared<ExpressionTree>();
        x->root = std::make_shared<ExpressionNode>();
        x->root->op = op;
        x->root->lhs = std::dynamic_pointer_cast<ExpressionTree>(variables.at(0))->root;
        x->root->rhs = std::dynamic_pointer_cast<ExpressionTree>(variables.at(2))->root;
        if(variables.size()==5) {
            x->root->rhs->tainted = true;
            x->root->rhs->counter = &counter;
            x->root->K = K;
        }
        return x;
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
    auto G = std::make_shared<StatefulLRParserBuilder>();
    std::string scope;
    int K,N;
    std::cin >> K >> N;
    int counter=0;
    G->addRuleList(std::make_shared<LastProjection>(),"Start","Expression");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Add), "Expression","Expression","+","Term");
    G->addRuleList(parser::reducers::Identity::instance, "Expression","Term");
//    G->addRuleList(parser::reducers::Identity::instance, "SEP","#");
//    G->addRuleList(parser::reducers::Identity::instance, "SEP","SEP","#");
    G->addRuleList(std::make_shared<ExpressionTreeBuilderWithCounter>(ExpressionNode::Add,counter,K), "Expression","Expression","+","Term","SEP");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Mul), "Term","Term","*","Factor");
    G->addRuleList(parser::reducers::Identity::instance, "Term","Factor");
    G->addRuleList(std::make_shared<NumberToLeaf>(), "Factor","Number");
    G->addRuleList(std::make_shared<ParenthesisBuilder>(), "Factor","(","Expression",")");
    G->addRuleList(std::make_shared<LeadingNumberConcatenation>(), "Number","NonZero","NumberWithZero"); // Number -> Number Digit
    G->addRuleList(std::make_shared<LeadingNumberConcatenation>(), "Number","NonZero"); // Number -> Number Digit
    G->addRuleList(std::make_shared<LeadingNumberConcatenation>(), "Number","Zero"); // Number -> Number Digit
    G->addRuleList(parser::reducers::Identity::instance, "NumberWithZero","Digit");
    G->addRuleList(std::make_shared<NumberConcatenation>(), "NumberWithZero","NumberWithZero","Digit");

    for(int i=0;i<10;i++)
        G->addRuleList(std::make_shared<NumberGenerator>(i), "Digit",std::to_string(i)); // Digit -> i
    for(int i=1;i<10;i++)
        G->addRuleList(std::make_shared<NumberGenerator>(i), "NonZero",std::to_string(i)); // Digit -> i
    G->addRuleList(std::make_shared<NumberGenerator>(0), "Zero","0"); // Digit -> i
    G->build();
    //G->printTable(std::cout);
    std::string expression;
    std::cin >> expression;
    for(int i=0;i<N;i++) for(int j=i+1;j<=N;j++)
    {
        auto sub=expression.substr(i,j-i);
        auto result = G->evaluate(sub);
        if(result)
        {
            auto main = std::dynamic_pointer_cast<ExpressionTree>(result);
            auto value = main->evaluate();
            counter += value%K==0;
            //std::cout << sub << " = " << value%K << '\n';
        }
    }
    std::cout << counter << '\n';

}


//int -> @
//class -> °
// "class x { int r;};" -> "° x { @ r;};"
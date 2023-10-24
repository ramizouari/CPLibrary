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

template<typename T>
struct Type : public parser::Variable
{
    T value;
    explicit Type(T value):value(std::move(value)){}
    Type() = default;
};
using StringValue = Type<std::string>;
using IntValue = Type<std::int64_t>;

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
    std::int64_t evaluate(class  Context & context);
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
    std::vector<std::map<std::string, std::int64_t>> variables_stack;
    std::map<std::string,std::vector<std::int64_t>> variables_index;
    std::map<std::string, std::shared_ptr<FunctionDefinition>> functions;

    Context() : variables_stack(1) {}

    void push()
    {
        variables_stack.emplace_back();
    }
    void pop()
    {
        for(auto &x: variables_stack.back())
        {
            variables_index[x.first].pop_back();
            if(variables_index[x.first].empty())
                variables_index.erase(x.first);
        }
        variables_stack.pop_back();
    }

    std::int64_t& get_variable(const std::string &name)
    {
        return variables_stack[variables_index.at(name).back()].at(name);
    }

    std::int64_t& upsert_variable(const std::string &name, std::int64_t value=0)
    {
        if(variables_index.find(name)==variables_index.end())
            variables_index[name]=std::vector<std::int64_t>(1,variables_stack.size()-1);
        else if(variables_index.at(name).back()!=variables_stack.size()-1)
            variables_index.at(name).push_back(variables_stack.size()-1);
        variables_stack.back()[name]=value;
        return variables_stack.back().at(name);
    }

    std::int64_t& get_variable(const std::string &name, int index)
    {
        return variables_stack[variables_index.at(name).at(index)].at(name);
    }


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
        Pow,
        Concat,
        Number,
        Variable,
        FunCall
    };
    Operation op;
    std::shared_ptr<ExpressionNode> lhs,rhs;
    std::variant<std::string, std::int64_t,FunctionCall> value;
    std::int64_t evaluate(Context &context);
};

struct ExpressionTree : public parser::Variable
{
    std::shared_ptr<ExpressionNode> root;
    std::int64_t evaluate(Context &context)
    {
        return root->evaluate(context);
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

std::int64_t ExpressionNode::evaluate(Context &context)
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
        case Pow:
            return power(lhs->evaluate(context),rhs->evaluate(context));
        case Number:
            return std::get<std::int64_t>(value);
        case Variable:
            return context.get_variable(std::get<std::string>(value));
        case FunCall:
        {
//            Context new_context;
//            new_context.functions = context.functions;
//            auto fn_call = std::get<FunctionCall>(value);
//            auto fn = new_context.functions.at(fn_call.name);
//            for(int i=0;i<fn->formal_parameters->value.size();i++)
//                new_context.variables[fn->formal_parameters->value[i]] = fn_call.parameters[i]->evaluate(context);
//            return fn->evaluate(new_context);
                auto fn_call = std::get<FunctionCall>(value);
                auto fn = context.functions.at(fn_call.name);
                std::vector<std::int64_t> parameters;
                parameters.reserve(fn_call.parameters.size());
                for(auto &x:fn_call.parameters)
                    parameters.push_back(x->evaluate(context));
                context.push();
                for(int i=0;i<fn->formal_parameters->value.size();i++)
                    context.upsert_variable(fn->formal_parameters->value[i],parameters[i]);
                auto result = fn->evaluate(context);
                context.pop();
                return result;
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
    std::int64_t evaluate(Context &context)
    {
//        auto result=context.variables[lhs] = rhs->evaluate(context);
        auto result=context.upsert_variable(lhs) = rhs->evaluate(context);
        if(lhs==OUTPUT_SYMBOL)
            std::cout << result << '\n';
        return result;
    }
};


struct ExecutionGraph : public parser::Variable
{
    std::vector<std::shared_ptr<Statement>> statements;
};

std::int64_t FunctionDefinition::evaluate(Context &context)
{
    std::int64_t result;
    for(auto statement : graph->statements)
        result = statement->evaluate(context);
    return result;
}

using Number = Type<std::int64_t>;

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
    std::shared_ptr<ExecutionGraph> main;
    std::shared_ptr<Context> context;
    std::int64_t evaluate()
    {
        std::int64_t result;
        for(auto statement : main->statements)
            result = statement->evaluate(*context);
        return result;
    }
};

struct PrepareExecution : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<MainProgram>();
        x->context = std::make_shared<Context>();
        if(v.size()==3)
        {
            x->main = std::dynamic_pointer_cast<ExecutionGraph>(v.at(2));
            x->context->functions = std::move(std::dynamic_pointer_cast<FunctionsBlock>(v.at(0))->functions);
        }
        else x->main = std::dynamic_pointer_cast<ExecutionGraph>(v.at(0));
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

struct OutputBuilder : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> &v) const override
    {
        auto x= std::make_shared<Statement>();
        x->lhs=OUTPUT_SYMBOL;
        x->rhs=std::dynamic_pointer_cast<ExpressionTree>(v.at(1));
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

class Encoder
{
    std::map<std::string,std::string> mapper;
    std::string alphabet;
    std::map<char,int> alphabet_map;
    std::string last_encoded;
public:
    explicit Encoder(std::string alphabet):alphabet(std::move(alphabet))
    {
        for(int i=0;i<this->alphabet.size();i++)
            alphabet_map[alphabet[i]]=i;
        mapper[""]="";
        mapper[OUTPUT_CODE] = OUTPUT_SYMBOL;
        mapper[RETURN_CODE] = RETURN_SYMBOL;
    }
    std::string encode(const std::string &x)
    {
        if(mapper.find(x)==mapper.end())
        {
            if(std::all_of(last_encoded.begin(),last_encoded.end(),[&alphabet=alphabet](auto x){return x==alphabet.back();}))
            {
                std::fill(last_encoded.begin(),last_encoded.end(),alphabet.front());
                last_encoded.push_back(alphabet.front());
                mapper[x]=last_encoded;
            }
            else
            {
                auto it=last_encoded.begin();
                while(*it==alphabet.back())
                {
                    *it=alphabet.front();
                    it++;
                }
                *it=alphabet[alphabet_map[*it]+1];
                mapper[x]=last_encoded;
            }
        }
        return mapper[x];
    }
};

bool isliteral(char x)
{
    return std::isalpha(x) || x == '_';
}

class Tokenizer
{
public:
    std::shared_ptr<parser::StatefulShiftReduceParser> builder;
    std::shared_ptr<Encoder> encoder;

    std::string tokenize(const std::string &source)
    {
        std::string result;
        result.reserve(source.size());
        int start,end;
        for(start=0;start<source.size() && std::isspace(source[start]);start++);
        for(end=source.size(); end > start && std::isspace(source[end-1]);end--);
        std::string word;
        for(int i=start;i< end; i++)
        {
            if(isliteral(source[i]))
                word.push_back(source[i]);
            else
            {
                result+=(encoder?encoder->encode(word):word);
                if(!std::isspace(source[i]) || source[i]=='\n')
                    result.push_back(source[i]);
                word.clear();
            }
        }
        return result;
    }
    auto evaluate(const std::string & source)
    {
        return builder->evaluate(tokenize(source));
    }
};

constexpr std::string_view alphabet="abcd";

/*
 * This is a simple interpreter for the VC language that supports the following features:
 * 1- Function definition
 * 2- Function call
 * 3- Variable assignment
 * 4- Output
 * 5- Return
 * 6- Arithmetic operations
 *
 * The VC language is composed of a potentially empty list of function definitions followed by a potentially empty list of output statements.
 * A function definition is composed of a function name, a list of formal parameters, and a function body.
 * A function body is composed of a potentially empty list of statements followed by a return statement.
 * A statement is either an assignment statement or an output statement.
 * */

int main(int argc, char** argv)
{
    using namespace parser;
    using namespace reducers;
    auto G = std::make_shared<StatefulLRParserBuilder>();
    std::string scope;

    G->addRuleList(std::make_shared<LastProjection>(),"Start","Program");
    G->addRuleList(std::make_shared<PrepareExecution>(),"Program","FunDefBlock","EOL", "Outputs");
    G->addRuleList(std::make_shared<PrepareExecution>(),"Program", "Outputs");
    G->addRuleList(std::make_shared<StatementConcat>(),"Outputs", "Outputs","EOL","Output");
    G->addRuleList(std::make_shared<ToListStatement>(),"Outputs", "Output");
    G->addRuleList(std::make_shared<OutputBuilder>(),"Output",OUTPUT_SYMBOL,"Expression");
    G->addRuleList(std::make_shared<FunctionsBlockConcat>(),"FunDefBlock","FunDefBlock","EOL","FunDef");
    G->addRuleList(std::make_shared<FunctionsBlockBuilder>(),"FunDefBlock","FunDef");
    G->addRuleList(std::make_shared<FnDefGenerator>(),"FunDef","Name","FunParams","FunBody");
    G->addRuleList(std::make_shared<Projection>(1),"FunParams","(","FunParamGroup",")","EOL");
    G->addRuleList(std::make_shared<DefaultGenerator<Type<std::vector<std::string>>>>(),"FunParams","(",")","EOL");
    G->addRuleList(std::make_shared<ListNames>(),"FunParamGroup","FunParamGroup","Comma","Name");
    G->addRuleList(std::make_shared<ToList<std::string>>(),"FunParamGroup","Name");
    G->addRuleList(reducers::Identity::instance,"EOL","EOL","\n");
    G->addRuleList(reducers::Identity::instance,"EOL","\n");
    G->addRuleList(std::make_shared<PrepareFunctionBody>(),"FunBody","StatementsBlock","EOL","Return");
    G->addRuleList(std::make_shared<PrepareFunctionBody>(),"FunBody","Return");
    G->addRuleList(std::make_shared<LastProjection>(),"Return",RETURN_SYMBOL,"Expression");
    G->addRuleList(std::make_shared<StatementConcat>(), "StatementsBlock","StatementsBlock","EOL","Statement");
    G->addRuleList(std::make_shared<ToListStatement>(), "StatementsBlock","Statement");
    G->addRuleList(parser::reducers::Identity::instance, "Comma",",");
    G->addRuleList(std::make_shared<StatementBuilder>(), "Statement","Name","=","Expression");
    G->addRuleList(reducers::Identity::instance, "Statement","Output");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Add), "Expression","Expression","+","Term");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Sub), "Expression","Expression","-","Term");
    G->addRuleList(parser::reducers::Identity::instance, "Expression","Term");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Mul), "Term","Term","*","Exponent");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Div), "Term","Term","/","Exponent");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Mod), "Term","Term","%","Exponent");
    G->addRuleList(parser::reducers::Identity::instance, "Term","Exponent");
    G->addRuleList(parser::reducers::Identity::instance, "Exponent","Factor");
    G->addRuleList(std::make_shared<ExpressionTreeBuilder>(ExpressionNode::Pow), "Exponent","Factor","^","Exponent");
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
    for(auto x:alphabet)
    {
        G->addRuleList(std::make_shared<CharacterGenerator>(x), "Name", std::string(1, x));
        G->addRuleList(std::make_shared<NameConcat>(x),"Name","Name",std::string(1, x));
    }
    G->build();
    G->printTable(std::cout);
    Tokenizer T;
    T.builder = G;
    T.encoder = std::make_shared<Encoder>(alphabet.data());
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
        main->evaluate();
    }
    else
    {
        std::cout << "error\n";
    }
}
//
// Created by ramizouari on 20/10/23.
//
#include "parser/StatefulParser.h"
#include "parser/LRParserBuilder.h"
#include <iostream>
#include <cmath>

namespace parser
{

}


struct Number : public parser::Variable
{
    std::int64_t v;
    Number(std::int64_t v):v(v){}
};

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

struct NumberAddition : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(std::dynamic_pointer_cast<Number>(variables.at(0))->v + std::dynamic_pointer_cast<Number>(variables.at(2))->v);
    }
};

struct NumberSubtraction: public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(std::dynamic_pointer_cast<Number>(variables.at(0))->v - std::dynamic_pointer_cast<Number>(variables.at(2))->v);
    }
};

struct NumberMultiplication: public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(std::dynamic_pointer_cast<Number>(variables.at(0))->v * std::dynamic_pointer_cast<Number>(variables.at(2))->v);
    }
};

struct NumberEuclideanDivision: public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(std::dynamic_pointer_cast<Number>(variables.at(0))->v / std::dynamic_pointer_cast<Number>(variables.at(2))->v);
    }
};
struct NumberModulo: public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(std::dynamic_pointer_cast<Number>(variables.at(0))->v % std::dynamic_pointer_cast<Number>(variables.at(2))->v);
    }
};

struct NumberConcatenation : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(10*std::dynamic_pointer_cast<Number>(variables.at(0))->v + std::dynamic_pointer_cast<Number>(variables.at(1))->v);
    }
};

struct Identity : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return variables.front();
    }
};

struct Inverse : public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(-std::dynamic_pointer_cast<Number>(variables.at(1))->v);
    }
};

struct Projection : public parser::VariableReducer
{
    int index;
    Projection(int index):index(index){}
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return variables.at(index);
    }
};

struct NumberExponentiation: public parser::VariableReducer
{
    [[nodiscard]] std::shared_ptr<parser::Variable> combine(const std::vector<std::shared_ptr<parser::Variable>> & variables) const override
    {
        return std::make_shared<Number>(std::pow(std::dynamic_pointer_cast<Number>(variables.at(0))->v, std::dynamic_pointer_cast<Number>(variables.at(2))->v));
    }
};

int main(int argc, char** argv)
{
    using namespace parser;
    StatefulLRParserBuilder LR;
    if(argc>=2)
    {
        auto grammar=std::string_view(argv[1]);
        if(grammar == "minimal")
        {

            LR.addRuleList(std::make_shared<Identity>(),"Start","E");
            LR.addRuleList(std::make_shared<Projection>(1),"E","(","E",")");
            LR.addRuleList(std::make_shared<NumberGenerator>(1),"E","1");
        }
        else if(grammar=="example")
        {
            LR.addRuleList(std::make_shared<Identity>(),"S","E");
            LR.addRuleList(std::make_shared<Identity>(),"E","T");
            LR.addRuleList(std::make_shared<Identity>(),"E","(","E",")");
            LR.addRuleList(std::make_shared<Identity>(),"T","n");
            LR.addRuleList(std::make_shared<Identity>(),"T","+","T");
            LR.addRuleList(std::make_shared<Identity>(),"T","T","+","n");

        }
        else
        {
            LR.addRuleList(std::make_shared<Identity>(),"Start","E");
            LR.addRuleList(std::make_shared<NumberAddition>(),"E","E","+","D");
            LR.addRuleList(std::make_shared<Identity>(),"E","D");
            LR.addRuleList(std::make_shared<Projection>(1),"D","(","E",")");
            LR.addRuleList(std::make_shared<NumberGenerator>(1),"D","1");
        }
    }
    else
    {
        LR.addRuleList(std::make_shared<Identity>(),"Start","E");
        LR.addRuleList(std::make_shared<NumberAddition>(),"E","E","+","M");
        LR.addRuleList(std::make_shared<NumberSubtraction>(),"E","E","-","M");
        LR.addRuleList(std::make_shared<Identity>(),"E","M");
        LR.addRuleList(std::make_shared<NumberMultiplication>(),"M","M","*","P");
        LR.addRuleList(std::make_shared<NumberEuclideanDivision>(),"M","M","/","P");
        LR.addRuleList(std::make_shared<NumberModulo>(),"M","M","%","P");
        LR.addRuleList(std::make_shared<Identity>(),"M","P");
        LR.addRuleList(std::make_shared<NumberExponentiation>(),"P","P","^","N");
        LR.addRuleList(std::make_shared<Identity>(),"P","N");

        LR.addRuleList(std::make_shared<Projection>(1),"N","(","B",")");
        LR.addRuleList(std::make_shared<Inverse>(),"B","-","E");
        LR.addRuleList(std::make_shared<Identity>(),"B","E");

        for(int i=0;i<10;i++)
        {
            LR.addRuleList(std::make_shared<NumberConcatenation>(), "N", "N","D");
            LR.addRuleList(std::make_shared<Identity>(), "N","D");
            LR.addRuleList(std::make_shared<NumberGenerator>(i), "D",std::to_string(i));
        }
    }
    LR.build();
    LR.printTable(std::cout);
    auto result=LR.evaluate("6+1+(2*31)-8^2");
    auto R=std::dynamic_pointer_cast<Number>(result);
    if(R)
        std::cout << R->v;
    else
    {
        std::cout << "Failed to calculate grammar";
    }
}
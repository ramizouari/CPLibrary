//
// Created by ramizouari on 20/10/23.
//
#include "parser/LRLogic.h"
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

int main()
{
    using namespace parser;
    LRLogic LR;
    LR.addFunctionalRuleList(std::make_shared<Identity>(),"S","Exp");
    LR.addFunctionalRuleList(std::make_shared<NumberExponentiation>(),"Exp","Exp","^","M");
    LR.addFunctionalRuleList(std::make_shared<Identity>(),"Exp","M");
    LR.addFunctionalRuleList(std::make_shared<NumberMultiplication>(),"M","M","*","E");
    LR.addFunctionalRuleList(std::make_shared<NumberEuclideanDivision>(),"M","M","/","E");
    LR.addFunctionalRuleList(std::make_shared<NumberModulo>(),"M","M","%","E");
    LR.addFunctionalRuleList(std::make_shared<Identity>(),"M","E");
    LR.addFunctionalRuleList(std::make_shared<NumberAddition>(),"E","E","+","N");
    LR.addFunctionalRuleList(std::make_shared<NumberSubtraction>(),"E","E","-","N");
    LR.addFunctionalRuleList(std::make_shared<Identity>(),"E","N");
    LR.addFunctionalRuleList(std::make_shared<Projection>(1),"N","(","B",")");
    LR.addFunctionalRuleList(std::make_shared<Inverse>(),"B","-","Exp");
    LR.addFunctionalRuleList(std::make_shared<Identity>(),"B","Exp");

    for(int i=0;i<10;i++)
    {
        LR.addFunctionalRuleList(std::make_shared<NumberGenerator>(i), "N", std::to_string(i));
        LR.addFunctionalRuleList(std::make_shared<NumberConcatenation>(), "N", "N",std::to_string(i));
    }
    LR.build();
    auto result=LR.calculate("5+((3-4)*0)+8+(5*2)+(((-3)*(-3))^2)");
    auto R=std::dynamic_pointer_cast<Number>(result);
    if(R)
        std::cout << R->v;
    else
    {
        std::cout << "Failed to calculate grammar";
    }
}
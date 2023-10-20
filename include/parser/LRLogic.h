//
// Created by ramizouari on 20/10/23.
//

#ifndef CPLIBRARY_LRLOGIC_H
#define CPLIBRARY_LRLOGIC_H
#include "LRParserBuilder.h"

namespace parser
{
    struct Variable
    {
        virtual ~Variable() = default;
    };

    struct VariableGenerator
    {
        std::shared_ptr<Variable> generate(Symbol &symbol) const;
    };

    struct VariableReducer
    {
        virtual ~VariableReducer() = default;
        virtual std::shared_ptr<Variable> combine(const std::vector<std::shared_ptr<Variable>>& symbols) const = 0;
    };

    struct NullReducer : public VariableReducer
    {
        std::shared_ptr<Variable> combine(const std::vector<std::shared_ptr<Variable>>& symbols) const
        {
            return nullptr;
        }
    };

    class LRLogic : public LRParserBuilder
    {
        std::vector<std::shared_ptr<VariableReducer>> reducers;
        inline static std::shared_ptr<VariableReducer> nullReducer = std::make_shared<NullReducer>();
        std::shared_ptr<Variable> result;
    public:
        using Action=LRFamily::Action;
        template<typename Left, typename ...Right>
        void addFunctionalRuleList(std::shared_ptr<VariableReducer> reducer, Left && left, Right && ... right)
        {
            this->LRParserBuilder::addRuleList(std::forward<Left>(left),std::forward<Right>(right)...);
            this->reducers.push_back(reducer);
        }

        template<typename Left, typename ...Right>
        void addRuleList(Left && left, Right&& ... right)
        {
            this->addFunctionalRuleList(nullReducer,std::forward<Left>(left),std::forward<Right>(right)...);
        }



        void addRule(std::shared_ptr<VariableReducer> reducer,const std::string &line)
        {
            LRParserBuilder::addRule(line);
            reducers.push_back(reducer);
        }

        void addRule(const std::string &line) override
        {
            this->addRule(nullReducer,line);
        }

        void addRuleList(std::shared_ptr<VariableReducer> reducer, const stringType &left,const std::vector<stringType> &right)
        {
            LRParserBuilder::addRuleList(left,right);
            reducers.push_back(reducer);
        }

        bool parse(const std::string &s) override;
        std::shared_ptr<Variable> calculate(const std::string &s);
    };
}

#endif //CPLIBRARY_LRLOGIC_H

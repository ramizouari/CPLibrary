//
// Created by ramizouari on 20/10/23.
//

#ifndef CPLIBRARY_STATEFULPARSER_H
#define CPLIBRARY_STATEFULPARSER_H
#include "Grammar.h"
#include "LRParserBuilder.h"

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
            inline static std::shared_ptr<VariableReducer> instance = std::make_shared<NullReducer>();
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

    class StatefulShiftReduceParser : virtual public StatefulStringParser<char>
    {
    public:
        using Action=LRFamily::Action;
        void printTable(std::ostream &H) const;
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
}

#endif //CPLIBRARY_STATEFULPARSER_H

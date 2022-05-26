#include <iostream>
#include <vector>
#include "parser/Grammar.h"
#include "parser/LRParser.h"

int main()
{
    parser::SLRParser parser;
    parser.addRuleList("E","E","+","M");
    parser.addRuleList("E","M");
    parser.addRuleList("M","M","*","T");
    parser.addRuleList("M","T");
    parser.addRuleList("T","(","E",")");
    parser.addRuleList("T","b","T");
    parser.addRuleList("T","a","T");
    parser.addRuleList("T","a");
    parser.addRuleList("T","b");

    parser.build();
    std::cout << parser.parse("(aaaa+bbababaab)+a+b*b*b*b*((a+b+a+(((b)*a)*bb)+abb))") << std::endl;
}
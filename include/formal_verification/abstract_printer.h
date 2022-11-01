//
// Created by ramizouari on 31/10/22.
//

#ifndef VFSF_ABSTRACT_PRINTER_H
#define VFSF_ABSTRACT_PRINTER_H

#include <ostream>
#include "dynamic_decision_tree.h"

namespace FormalSpecification::Printer {

    struct abstract_printer
    {
        std::ostream &ostream;
        abstract_printer(std::ostream &S);
    protected:
        virtual void print(const dynamic_decision_tree &T)=0;
        friend abstract_printer& operator<<(abstract_printer&P,const dynamic_decision_tree&T);
    };

    template<typename T>
    concept Printable= requires(std::ostream &H,const T&x)
    {
        {H << x} -> std::convertible_to<std::ostream&>;
    };


    struct mermaid_printer:public abstract_printer
    {
        using abstract_printer::abstract_printer;
        void print(const dynamic_decision_tree &T) override;
    };

    struct graphml_printer:public abstract_printer
    {
        using abstract_printer::abstract_printer;
        void print(const dynamic_decision_tree &T) override;
    };

    abstract_printer& operator<<(abstract_printer&P, const dynamic_decision_tree&T);

    template<Printable T>
    abstract_printer& operator<<(abstract_printer& P,const T&x)
    {
        P.ostream << x;
        return P;
    }
} // Printer

#endif //VFSF_ABSTRACT_PRINTER_H

#ifndef UTOPIA_FE_VAR_HPP
#define UTOPIA_FE_VAR_HPP

#include "utopia_Options.hpp"

#include <iostream>
#include <string>

namespace utopia {

    class FEVar : public Configurable, public Describable {
    public:
        void read(Input &in) override {
            if (!Options()
                     .add_option("fe_family", fe_family, "The finite element family from enum libMesh::FEFamily.")
                     .add_option("order", order, "The finite element order from enum libMesh::Order.")
                     .add_option("name", name, "Name of the variable")
                     .add_option("n_components", n_components, "number of vector n_component")
                     .parse(in)) {
                return;
            }
        }

        void describe(std::ostream &os) const override {
            os << "fe_family:\t" << fe_family;
            os << "order:\t" << order;
            os << "name:\t" << name;
            os << "n_components:\t" << n_components;
        }

        std::string fe_family{"LAGRANGE"};
        std::string order{"FIRST"};
        std::string name{"var"};
        int n_components{1};
    };

}  // namespace utopia

#endif  // UTOPIA_FE_VAR_HPP
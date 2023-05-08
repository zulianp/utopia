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
        std::string name{"fevar"};
        int n_components{1};
    };

    class FEVariables : public Configurable, public Describable {
    public:
        using ConstIterator = std::vector<FEVar>::const_iterator;

        void read(Input &in) override {
            int n_var{1};

            if (!Options().add_option("n_var", n_var, "number of vector n_component").parse(in)) {
                return;
            }

            in.get("variables", [this](Input &in) {
                in.get_all([&](Input &in) {
                    FEVar v;
                    v.read(in);
                    variables_.push_back(v);
                });
            });

            if (variables_.empty()) {
                FEVar v;
                // At least one variable
                v.n_components = std::max(1, n_var);
                variables_.push_back(v);
            }
        }

        int count_variables() {
            int counted_vars = 0;

            for (auto &v : variables_) {
                counted_vars += v.n_components;
            }

            return counted_vars;
        }

        ConstIterator begin() const { return variables_.begin(); }
        ConstIterator end() const { return variables_.end(); }

        FEVar &operator[](const std::size_t i) { return variables_[i]; }
        const FEVar &operator[](const std::size_t i) const { return variables_[i]; }
        bool empty() const { return variables_.empty(); }

        void describe(std::ostream &os) const override {
            for (auto &v : variables_) {
                v.describe(os);
            }
        }

    private:
        std::vector<FEVar> variables_;
    };

}  // namespace utopia

#endif  // UTOPIA_FE_VAR_HPP
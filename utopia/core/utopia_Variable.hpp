#ifndef UTOPIA_VARIABLE_HPP
#define UTOPIA_VARIABLE_HPP

#include <memory>
#include "utopia_Expression.hpp"
#include "utopia_TreeNavigator.hpp"

namespace utopia {

    template <class Expr, int Number>
    class Variable : public Expression<Variable<Expr, Number> > {
    public:
        static const int Order = Expr::Order;

        using Scalar = typename Expr::Scalar;

        inline Expr &expr() {
            assert(expr_);
            return *expr_;
        }

        inline const Expr &expr() const {
            assert(expr_);
            return *expr_;
        }

        inline static constexpr int get_number() { return Number; }

        void set(const std::shared_ptr<Expr> &expr) { expr_ = expr; }

        Variable(const std::shared_ptr<Expr> &expr) : expr_(expr) {}

        std::string get_class() const override { return "Variable<" + expr_->get_class() + ">"; }

    private:
        std::shared_ptr<Expr> expr_;
    };

    template <class Expr, int Number>
    class SearchVariableAction {
    public:
        template <class Node>
        inline void pre_order_visit(const Node &) {}

        template <class Node>
        inline void post_order_visit(const Node &) {}

        template <class Node>
        inline void in_order_visit(const Node &) {}

        void pre_order_visit(const Variable<Expr, Number> &var) { var_ = &var; }

        SearchVariableAction() : var_(nullptr) {}

        const Variable<Expr, Number> *get_var() {
            assert(var_);
            return var_;
        }

    public:
        const Variable<Expr, Number> *var_;
    };

    template <class VarExpr, int Number, class Expr>
    Variable<VarExpr, Number> &find_variable(Expression<Expr> &expr) {
        // TODO make if more efficient by using a static search and pruning
        SearchVariableAction<VarExpr, Number> action;
        auto nav = make_nav(action);
        nav.visit(expr.derived());
        // this is not so bad because of the signature.
        return *const_cast<Variable<VarExpr, Number> *>(action.get_var());
    }

    template <class Expr, int Number>
    class Traits<Variable<Expr, Number> > : public Traits<Expr> {};

    template <class Expr, int Number>
    auto size(const Variable<Expr, Number> &expr) -> decltype(size(expr.expr())) {
        return size(expr.expr());
    }

    template <class Expr, int Number, class Traits, int Backend>
    class Eval<Variable<Expr, Number>, Traits, Backend> {
    public:
        inline static auto apply(const Variable<Expr, Number> &expr) -> decltype(expr.expr().derived()) {
            return expr.expr().derived();
        }
    };
}  // namespace utopia

#endif  // UTOPIA_VARIABLE_HPP

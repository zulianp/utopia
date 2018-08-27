//
// Created by Patrick Zulian on 31/08/16.
//

#ifndef UTOPIA_UTOPIA_PRECONDITIONER_HPP
#define UTOPIA_UTOPIA_PRECONDITIONER_HPP

#include "utopia_Parameters.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Expression.hpp"

#include <memory>

namespace utopia {
    template<class Vector>
    class Preconditioner {
    public:
        virtual ~Preconditioner() {}
        virtual bool apply(const Vector &rhs, Vector &sol) = 0;
        virtual void set_parameters(const Parameters)
        {}
    };



    template<class Expr, class Vector>
    class ExprPreconditioner : public Preconditioner<Vector> {
    public:
        bool apply(const Vector &rhs, Vector &sol) override
        {
            sol = expr_ * rhs;
            std::cout<< "expression thing ... \n";
            return true;
        }

        ExprPreconditioner(const Expr &expr)
        : expr_(expr)
        {}

    private:
        UTOPIA_STORE_CONST(Expr) expr_;
    };

    template<class Vector, class Derived>
    std::shared_ptr<ExprPreconditioner<Derived, Vector> > make_preconditioner(const Expression<Derived> &expr)
    {
        return std::make_shared<ExprPreconditioner<Derived, Vector> >(expr.derived());
    }

    template<class Matrix, class Vector>
    class DelegatePreconditioner : public Preconditioner<Vector> {
    public:

        bool apply(const Vector &/*rhs*/, Vector &/*sol*/) override
        {
            std::cerr<< "[Warning] DelegatePreconditioner::apply doing nothing ... \n";
            return true;
        }

        void update(const std::shared_ptr<const Matrix> &op)
        {
            // disp(*op);
            // std::cout<< "[Message] DelegatePreconditioner::update ... \n";
            op_ = op;
        }

        const std::shared_ptr<const Matrix> &get_matrix() const
        {
            // std::cout<< "[Message] DelegatePreconditioner::get_matrix ... \n";
            return op_;
        }

    private:
        std::shared_ptr<const Matrix> op_;
    };
}

#endif //UTOPIA_UTOPIA_PRECONDITIONER_HPP

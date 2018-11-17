#ifndef UTOPIA_UTOPIA_PRECONDITIONER_HPP
#define UTOPIA_UTOPIA_PRECONDITIONER_HPP

#include "utopia_Parameters.hpp"
#include "utopia_StoreAs.hpp"
#include "utopia_Expression.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Traits.hpp"

#include <memory>
#include <cassert>

#define UTOPIA_W_VECTOR(Tensor) utopia::Wrapper<typename utopia::Traits<Tensor>::Vector, 1>

namespace utopia 
{
    template<class Vector>
    class Operator {
    public:
        virtual ~Operator() {}
        virtual bool apply(const Vector &rhs, Vector &sol) const = 0;
    };

    template<class Matrix, class Vector>
    class MatrixOperator final : public Operator<Vector> {
    public:
        MatrixOperator(const std::shared_ptr<const Matrix> &mat)
        : mat_(mat)
        {}

        bool apply(const Vector &rhs, Vector &ret) const override
        {
            ret = (*mat_) * rhs;
            return true;
        }

       const std::shared_ptr<const Matrix> &get_matrix() const
       {
           return mat_;
       }

    private:
        std::shared_ptr<const Matrix> mat_;
    };

    template<class Matrix>
    std::unique_ptr<
        MatrixOperator<Matrix, UTOPIA_W_VECTOR(Matrix)>
    > op(const std::shared_ptr<const Matrix> &mat)
    {
        return utopia::make_unique< MatrixOperator<Matrix, UTOPIA_W_VECTOR(Matrix)> >(mat);
    }

    template<class Matrix>
    std::unique_ptr<
        MatrixOperator<Matrix, UTOPIA_W_VECTOR(Matrix)>
    > op(const std::shared_ptr<Matrix> &mat)
    {
        return utopia::make_unique< MatrixOperator<Matrix, UTOPIA_W_VECTOR(Matrix)> >(mat);
    }


    template<class Vector>
    class Preconditioner /*: public Operator<Vector>*/ {
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
            op_ = op;
        }

        const std::shared_ptr<const Matrix> &get_matrix() const
        {
            return op_;
        }

    private:
        std::shared_ptr<const Matrix> op_;
    };
    
}

#endif //UTOPIA_UTOPIA_PRECONDITIONER_HPP

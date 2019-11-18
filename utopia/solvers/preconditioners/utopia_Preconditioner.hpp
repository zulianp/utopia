#ifndef UTOPIA_UTOPIA_PRECONDITIONER_HPP
#define UTOPIA_UTOPIA_PRECONDITIONER_HPP

#include "utopia_StoreAs.hpp"
#include "utopia_Expression.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_Traits.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_Input.hpp"
#include "utopia_Operator.hpp"
#include "utopia_Communicator.hpp"

#include <memory>
#include <cassert>

// #define UTOPIA_W_VECTOR(Tensor) utopia::Wrapper<typename utopia::Traits<Tensor>::Vector, 1>

namespace utopia {

    // template<class Matrix, class Vector>
    // class MatrixOperator final : public Operator<Vector> {
    // public:
    //     using Communicator = typename Traits<Vector>::Communicator;

    //     MatrixOperator(const std::shared_ptr<Matrix> &mat)
    //     : mat_(mat)
    //     {}

    //     bool apply(const Vector &rhs, Vector &ret) const override
    //     {
    //         ret = (*mat_) * rhs;
    //         return true;
    //     }

    //    // const std::shared_ptr<const Matrix> &get_matrix() const
    //    // {
    //    //     return mat_;
    //    // }

    //    inline Size size() const override
    //    {
    //         assert(mat_);
    //         return mat_->size();
    //    }

    //    inline Size local_size() const
    //    {
    //         assert(mat_);
    //         return mat_->local_size();
    //    }

    //     Communicator &comm()
    //    {
    //         assert(mat_);
    //         return mat_->comm();
    //    }

    //     const Communicator &comm() const
    //    {
    //         assert(mat_);
    //         return mat_->comm();
    //    }


    // private:
    //     std::shared_ptr<Matrix> mat_;
    // };


    template<class Vector, class Fun>
    class LambdaOperator final : public Operator<Vector> {
    public:
        using Communicator = typename Traits<Vector>::Communicator;

        LambdaOperator(
            Communicator &comm,
            const Size &size,
            const Size &local_size,
            Fun fun) :
        comm_(comm),
        size_(size),
        local_size_(local_size),
        fun_(fun)
        {}

        inline bool apply(const Vector &rhs, Vector &sol) const override
        {
            return fun_(rhs, sol);
        }

        inline Size size() const override
        {
            return size_;
        }

        inline Size local_size() const override
        {
            return local_size_;
        }

        inline Communicator &comm() override
        {
            return comm_;
        }

        inline const Communicator &comm() const override
        {
            return comm_;
        }

    private:
        Communicator &comm_;
        Size size_, local_size_;
        Fun fun_;
    };

    template<typename Vector>
    std::unique_ptr< LambdaOperator<Vector, std::function<bool(const Vector &, Vector &)>> > op(
        typename Traits<Vector>::Communicator &comm,
        const Size &size,
        const Size &local_size,
        std::function<bool(const Vector &, Vector &)> f)
    {
        return utopia::make_unique< LambdaOperator<Vector, std::function<bool(const Vector &, Vector &)>> >(
            comm,
            size,
            local_size,
            f
        );
    }

    // template<class Matrix>
    // std::unique_ptr< MatrixOperator<Matrix, UTOPIA_VECTOR(Matrix)> > op(const std::shared_ptr<const Matrix> &mat)
    // {
    //     return utopia::make_unique< MatrixOperator<Matrix, UTOPIA_VECTOR(Matrix)> >(mat);
    // }

    // template<class Matrix>
    // std::unique_ptr< MatrixOperator<Matrix, UTOPIA_VECTOR(Matrix)> > op_ref(const Matrix &mat)
    // {
    //     return utopia::make_unique< MatrixOperator<Matrix, UTOPIA_VECTOR(Matrix)> >(make_ref(mat));
    // }


    template<class Vector>
    class Preconditioner : public virtual Configurable, public virtual Clonable {
    public:
        virtual ~Preconditioner() {}
        virtual bool apply(const Vector &rhs, Vector &sol) = 0;

        virtual void read(Input &/*in*/) override
        {

        }

        virtual void print_usage(std::ostream & /*os*/) const override
        {

        }

        // virtual void update(const Operator<Vector> &) {}

        virtual Preconditioner * clone() const override = 0;
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

        virtual ExprPreconditioner * clone() const override
        {
            return new ExprPreconditioner(*this);
        }

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

        virtual DelegatePreconditioner * clone() const override
        {
            return new DelegatePreconditioner(*this);
        }

    private:
        std::shared_ptr<const Matrix> op_;
    };

}

#endif //UTOPIA_UTOPIA_PRECONDITIONER_HPP

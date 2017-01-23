//
// Created by Patrick Zulian on 31/08/16.
//

#ifndef UTOPIA_UTOPIA_PRECONDITIONER_HPP
#define UTOPIA_UTOPIA_PRECONDITIONER_HPP

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
            std::cout<< "expresion thing ... \n"; 
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
}

#endif //UTOPIA_UTOPIA_PRECONDITIONER_HPP

#ifndef UTOPIA_QUADRATIC_FUNCTION_HPP
#define UTOPIA_QUADRATIC_FUNCTION_HPP

#include "utopia_Function.hpp"

namespace utopia {
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class QuadraticFunction final : public Function<Matrix, Vector, Backend> {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        QuadraticFunction(const std::shared_ptr<Matrix> &H, const std::shared_ptr<Vector> &rhs)
        : rhs_(rhs)
        {
            assert(H);
            assert(rhs);
            
            this->data()->H = H;
        }

        virtual bool initialize_hessian(Matrix &H, Matrix &/*H_pre*/) const override
        {
            H = *this->data()->H;
            return true;
        }

        virtual ~QuadraticFunction() { }

        virtual bool value(const Vector &x, Scalar &value) const override
        {
            value = 0.5 * dot(x, (*this->data()->H) * x) - dot(x,  *rhs_);
            return true;
        }

        virtual bool gradient(const Vector &x, Vector &result) const override
        {
            result = (*this->data()->H) * x - (*rhs_);
            return true;
        }

        virtual bool hessian(const Vector &x, Matrix &H) const override
        {
            H = *this->data()->H;
            return true;
        }

        virtual bool hessian(const Vector &x, Matrix &result, Matrix &prec) const override
        {
            return false;
        }

        virtual bool has_preconditioner() const override
        {
            return false;
        }

        virtual bool update(const Vector &x) override
        {
            return true;
        }

    private:
        std::shared_ptr<Vector> rhs_;
    };
}

#endif //UTOPIA_QUADRATIC_FUNCTION_HPP

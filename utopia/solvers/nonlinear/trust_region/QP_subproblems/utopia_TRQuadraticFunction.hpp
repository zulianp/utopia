#ifndef UTOPIA_TR_QUADRATIC_FUNCTION_HPP
#define UTOPIA_TR_QUADRATIC_FUNCTION_HPP

#include "utopia_Function.hpp"

namespace utopia {

    // this function differs from QuadraticFunction in 2 things: - different signs for rhs term
    //                                                           - supports const
    // TODO:: merge with original QPFunction...                                                           
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class TRQuadraticFunction final : public Function<Matrix, Vector, Backend> {
    public:
        DEF_UTOPIA_SCALAR(Matrix)

        TRQuadraticFunction(const std::shared_ptr<const Matrix> &H, const std::shared_ptr<const Vector> &rhs)
        : rhs_(rhs), H_(H)
        {

        }
        
        ~TRQuadraticFunction() { }

        bool value(const Vector &x, Scalar &value) const override
        {
            value = 0.5 * dot(x, *H_ * x); 
            value += dot(x,  *rhs_);
            return true;
        }

        bool gradient(const Vector &x, Vector &result) const override
        {
            result = *H_ * x + (*rhs_);
            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            H = *H_;
            return true;
        }

        bool hessian(const Vector & /*x*/, Matrix &/*result*/, Matrix &/*prec*/) const override
        {
            return false;
        }

        bool has_preconditioner() const override
        {
            return false;
        }

        bool update(const Vector &/*x*/) override
        {
            return true;
        }

        bool initialize_hessian(Matrix &H, Matrix &/*H_pre*/) const override
        {
            H = *H_;
            return true;
        }


    private:
        std::shared_ptr<const Vector> rhs_;
        std::shared_ptr<const Matrix> H_;
    };
}

#endif //UTOPIA_QUADRATIC_FUNCTION_HPP

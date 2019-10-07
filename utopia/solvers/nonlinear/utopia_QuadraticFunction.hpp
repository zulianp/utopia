#ifndef UTOPIA_QUADRATIC_FUNCTION_HPP
#define UTOPIA_QUADRATIC_FUNCTION_HPP

#include "utopia_Function.hpp"

namespace utopia {
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class QuadraticFunction final : public Function<Matrix, Vector, Backend> {
    public:
        DEF_UTOPIA_SCALAR(Matrix);

        QuadraticFunction(const std::shared_ptr<Matrix> &H, const std::shared_ptr<Vector> &rhs)
        : rhs_(rhs), H_(H)
        {
            assert(H);
            assert(rhs);

            // this->data()->H = H;
        }

        bool initialize_hessian(Matrix &H, Matrix &/*H_pre*/) const override
        {
            // H = *this->data()->H;
            H = *H_;
            return true;
        }

        bool value(const Vector &x, Scalar &value) const override
        {
            value = 0.5 * dot(x, (*H_) * x) - dot(x, *rhs_);
            return true;
        }

        bool gradient(const Vector &x, Vector &result) const override
        {
            result = (*H_) * x - (*rhs_);
            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            UTOPIA_UNUSED(x);
            // H = *this->data()->H;
            H = *H_;
            return true;
        }

        bool hessian(const Vector &x, Matrix &result, Matrix &prec) const override
        {
            UTOPIA_UNUSED(x);
            UTOPIA_UNUSED(result);
            UTOPIA_UNUSED(prec);
            return false;
        }

        bool has_preconditioner() const override
        {
            return false;
        }

        bool update(const Vector &x) override
        {
            UTOPIA_UNUSED(x);
            return true;
        }

    private:
        std::shared_ptr<Vector> rhs_;
        std::shared_ptr<Matrix> H_;
    };






    template<class Matrix, class Vector>
    class QuadraticExtendedFunction final : public ExtendedFunction<Matrix, Vector> {
    public:
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef UTOPIA_SCALAR(Vector) Scalar;

        QuadraticExtendedFunction(  const Matrix &H, 
                                    const Vector &rhs, 
                                    const Vector & x_0, 
                                    const Vector & bc_marker, 
                                    const Vector & rhs_non_eq): 
                                    ExtendedFunction<Matrix, Vector>(x_0, bc_marker, rhs_non_eq),  
                                    rhs_(rhs), 
                                    H_(H)
        {
            // this->data()->H = H;
        }

        bool initialize_hessian(Matrix &H, Matrix &/*H_pre*/) const override
        {
            // H = *this->data()->H;
            H = H_;
            return true;
        }

        bool value(const Vector &x, Scalar &value) const override
        {
            value = 0.5 * dot(x, H_ * x) + dot(x, rhs_);
            return true;
        }

        bool gradient_no_rhs(const Vector &x, Vector &result) const override
        {
            result = H_ * x+ rhs_;

            {
                Vector bc_flg = this->get_eq_constrains_flg(); 

                Read<Vector>  r_ub(bc_flg); 
                Write<Vector> wv(result);

                Range r = range(result);

                for(SizeType i = r.begin(); i != r.end(); ++i)
                {
                    if(bc_flg.get(i) > 0.0){
                        result.set(i, 0.0); 
                    }
                }                       
            }


            return true;
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            UTOPIA_UNUSED(x);
            
            H = H_;
            
            const std::vector<SizeType> & index = this->get_indices_related_to_BC(); 
            set_zero_rows(H, index, 1.);

            return true;
        }

        bool hessian(const Vector &x, Matrix &result, Matrix &prec) const override
        {
            UTOPIA_UNUSED(x);
            UTOPIA_UNUSED(result);
            UTOPIA_UNUSED(prec);
            return false;
        }

        bool has_preconditioner() const override
        {
            return false;
        }

        bool update(const Vector &x) override
        {
            UTOPIA_UNUSED(x);
            return true;
        }

    private:
        Vector rhs_;
        Matrix H_;

    };

}

#endif //UTOPIA_QUADRATIC_FUNCTION_HPP

#ifndef UTOPIA_BRATU_3D_MAT_BASED_HPP
#define UTOPIA_BRATU_3D_MAT_BASED_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"


namespace utopia
{
    template<class Matrix, class Vector>
    class Bratu3D : public ExtendedFunction<Matrix, Vector> 
    {
    public:
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using Scalar = typename utopia::Traits<Vector>::Scalar;

        Bratu3D(    const Matrix & H, 
                    const Vector & x_0, 
                    const Vector & bc_marker, 
                    const Vector & rhs_non_eq, 
                    const Scalar & lambda, 
                    const Scalar & dx, 
                    const Scalar & dy,
                    const Scalar & dz): 
                    ExtendedFunction<Matrix, Vector>(x_0, bc_marker, rhs_non_eq),  
                    H_(H), 
                    lambda_(lambda), 
                    dx_(dx), 
                    dy_(dy), 
                    dz_(dz)
        {
            dxyzLambda_ = dx_ * dy_ * dz_ * lambda_; 
        }

        bool initialize_hessian(Matrix &H, Matrix &/*H_pre*/) const override
        {
            H = H_;
            return true;
        }

        bool value(const Vector &x, Scalar &value) const override
        {
            value = 0.5 * dot(x, (H_) * x) - (dxyzLambda_ * sum(exp(x))); 
            return true;
        }

        bool gradient(const Vector &x, Vector &result) const override
        {
            result = (H_ * x) - (dxyzLambda_ * exp(x)); 

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
            
            Matrix D = diag(exp(x)); 
            H = H_ - (dxyzLambda_* D); 

            
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
        Scalar lambda_; 
        Scalar dx_; 
        Scalar dy_; 
        Scalar dz_; 
        Scalar dxyzLambda_; 
    };
}


#endif //UTOPIA_BRATU_3D_MAT_BASED_HPP
#ifndef UTOPIA_EXTENDED_FUNCTION_HPP
#define UTOPIA_EXTENDED_FUNCTION_HPP

#include "utopia_Base.hpp"

namespace utopia
{
    /**
     * @brief      Class for Nonlinear Function, all application context needed by solver is usually provided inside of this functions.
     *             In optimization settings, user needs to supply value(energy), gradient, hessian.
     *             Difference, between Function and ExtendedFunction is that here, we make use of additional informations to improve convergence
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector>
    class ExtendedFunction : public Function<Matrix, Vector>
    {
    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        virtual ~ExtendedFunction() { }

        ExtendedFunction() {}

        ExtendedFunction(const Vector & x_init, const Vector & bc_marker, const Vector & rhs) :
                _x_eq_values(x_init),
                _rhs(rhs),
                _eq_constrains_flg(bc_marker)
        {
            Vector ones = local_values(local_size(_eq_constrains_flg).get(0), 1.0); 
            _eq_constraints_mask_matrix_ = diag(ones - _eq_constrains_flg); 

            {
                Read<Vector> r(_eq_constrains_flg);

                Range range_w = range(_eq_constrains_flg);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++)
                {
                    if(_eq_constrains_flg.get(i) == 1)
                        indices_eq_constraints_.push_back(i);
                }
            }              
        }

        virtual bool value(const Vector &/*point*/, Scalar &/*value*/) const override = 0;


        bool gradient(const Vector & x, Vector & gradient) const override final
        {
            this->gradient_no_rhs(x, gradient);

            if(size(gradient) == size(this->_rhs)) {
                gradient = gradient - this->_rhs;
            }

            return true;
        }


        virtual bool gradient_no_rhs(const Vector &/*point*/, Vector &/*result*/) const = 0;


        virtual bool hessian(const Vector &x, Matrix &H) const override = 0;
        virtual bool hessian(const Vector &/*point*/, Matrix &/*result*/, Matrix &/*preconditioner*/) const  override
        {
            return false;
        }

        virtual bool has_preconditioner() const override
        {
            return false;
        }

        virtual bool update(const Vector &/*point*/)  override
        {
            return true;
        }

        virtual bool set_rhs(const Vector & rhs)
        {
            _rhs = rhs;
            return true;
        }

        virtual bool reset_rhs()
        {
            _rhs = local_zeros(local_size(_rhs));
            return true;
        }


        virtual bool get_rhs( Vector & rhs) const
        {
            rhs = _rhs;
            return true;
        }

        virtual bool has_rhs() const
        {
            return !empty(_rhs);
        }

        virtual bool get_eq_constrains_values(Vector & x) const 
        {
            x = _x_eq_values;
            return true;
        }

        virtual bool get_eq_constrains_flg(Vector & x) const
        {
            x = _eq_constrains_flg;
            return true;
        }

        inline const std::vector<SizeType> & get_indices_related_to_BC() const
        {
            return indices_eq_constraints_; 
        }


        const Vector &get_eq_constrains_flg() const
        {
            return _eq_constrains_flg;
        }

        virtual bool set_equality_constrains(const Vector &eq_constrains_flg, const Vector &x_in)
        {
            _x_eq_values        =  x_in;
            _eq_constrains_flg  = eq_constrains_flg;

            Vector ones = local_values(local_size(_eq_constrains_flg).get(0), 1.0); 
            _eq_constraints_mask_matrix_ = diag(ones - _eq_constrains_flg); 

            {
                Read<Vector> r(_eq_constrains_flg);

                Range range_w = range(_eq_constrains_flg);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++)
                {
                    if(_eq_constrains_flg.get(i) == 1)
                        indices_eq_constraints_.push_back(i);
                }
            }            

            return true;
        }


        virtual bool zero_contribution_to_equality_constrains(Vector & x) const
        {
            x = _eq_constraints_mask_matrix_ * x; 
            return true;
        }




     protected:
        Vector _x_eq_values;
        Vector _rhs;
        Vector _eq_constrains_flg;
        
        Matrix _eq_constraints_mask_matrix_; 
        std::vector<SizeType> indices_eq_constraints_; 


    };
}
#endif //UTOPIA_EXTENDED_FUNCTION_HPP

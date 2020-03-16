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

        ExtendedFunction(const Vector & x_init, const Vector & bc_marker)
        {
            this->set_equality_constrains(bc_marker, x_init);               
        }

        virtual bool value(const Vector &/*point*/, Scalar &/*value*/) const override = 0;

        // Copy of vec... 
        Vector initial_guess() const
        {   
            return _x_eq_values; 
        }

        virtual SizeType loc_size() const 
        {
            return local_size(_x_eq_values).get(0); 
        }

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


        Vector &get_eq_constrains_flg() 
        {
            return _eq_constrains_flg;
        }

        Vector &get_eq_constrains_values()  
        {
            return  _x_eq_values;
        }        

        virtual bool set_equality_constrains(const Vector &eq_constrains_flg, const Vector &x_in)
        {
            _x_eq_values        =  x_in;
            _eq_constrains_flg  = eq_constrains_flg;

            this->init_constraint_indices(); 

            return true;
        }


        bool init_constraint_indices()
        {
            indices_eq_constraints_.clear();

            {
                Read<Vector> r(_eq_constrains_flg);

                Range range_w = range(_eq_constrains_flg);
                for (SizeType i = range_w.begin(); i != range_w.end(); i++)
                {
                    if(_eq_constrains_flg.get(i) == 1){
                        indices_eq_constraints_.push_back(i);
                    }
                }
            }      

            return true;
        }



        virtual bool zero_contribution_to_equality_constrains(Vector & x) const
        {
            UTOPIA_NO_ALLOC_BEGIN("RMTR::zero_contribution_to_equality_constrains");
            // x = _eq_constraints_mask_matrix_ * x; 

            {
                auto d_flg     = const_device_view(_eq_constrains_flg);

                parallel_transform(
                          x,
                          UTOPIA_LAMBDA(const SizeType &i, const Scalar &xi) -> Scalar {
                            Scalar flg = d_flg.get(i);
                            
                            // TODO:: use abs with eps tolerance
                            if(flg==1.0){
                              return 0.0;
                            }
                            else
                              return xi; 
                      });
            }


            UTOPIA_NO_ALLOC_END();
            return true;
        }


     protected:
        Vector _x_eq_values;
        Vector _eq_constrains_flg;
        
        std::vector<SizeType> indices_eq_constraints_; 
    };

    template<class Matrix, class Vector>
    class Function_rhs final: public ExtendedFunction<Matrix, Vector>
    {
    public:
        DEF_UTOPIA_SCALAR(Matrix);
        typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;

        typedef utopia::ExtendedFunction<Matrix, Vector>    Fun;

        Function_rhs(const std::shared_ptr<Fun> & fun): fun_(fun)
        {

        }

        bool value(const Vector &x, Scalar &value) const override
        {
            fun_->value(x, value); 
            return true; 
        }

        bool gradient(const Vector & x, Vector &g) const override
        {
            fun_->gradient(x, g); 

            if(size(g) == size(this->rhs_)) {
                g = g - this->rhs_;
            }            

            return true; 
        }

        bool hessian(const Vector &x, Matrix &H) const override
        {
            fun_->hessian(x, H); 
            return true; 
        }

        virtual bool set_rhs(const Vector & rhs)
        {
            rhs_ = rhs;
            return true;
        }

        virtual bool reset_rhs()
        {
            if(!empty(rhs_)){
                rhs_.set(0.0);
            }
            else{
                utopia_error("error in reset rhs... \n"); 
            }
            return true;
        }


        virtual bool get_rhs( Vector & rhs) const
        {
            rhs = rhs_;
            return true;
        }

        virtual bool has_rhs() const
        {
            return !empty(rhs_);
        }        

        private:
            std::shared_ptr<Fun> fun_;   
            Vector rhs_;
    };


}
#endif //UTOPIA_EXTENDED_FUNCTION_HPP

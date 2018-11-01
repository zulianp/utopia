#ifndef UTOPIA_HESSIAN_BFGS_HPP
#define UTOPIA_HESSIAN_BFGS_HPP

#include "utopia_Core.hpp"


namespace utopia
{

    template<class Matrix, class Vector>
    class BFGS final: public HessianApproximation<Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        static_assert(!utopia::is_sparse<Matrix>::value, "utopia::BFGS does not support sparse matrices.");

        public:

            BFGS(): HessianApproximation<Vector>(), update_hessian_(false)
            {

            }

            inline BFGS<Matrix, Vector> * clone() const override
            {
                return new BFGS<Matrix, Vector>(*this);
            }

            virtual void initialize(const SizeType & n) override
            {
                if(update_hessian_)
                {
                    H_prev_ = local_identity(n, n);
                }

                H_prev_inv_ = local_identity(n, n);

                s_ = local_zeros(n); 
                y_ = local_zeros(n); 

                this->initialized(true); 
            }                

            virtual bool update(const Vector & s_in, const Vector & y_in ) override
            {
                if(!this->initialized())
                {
                    utopia_error("BFGS::update: Initialization needs to be done before updating. \n"); 
                    return false; 
                }


                s_ = s_in; 
                y_ = y_in;
                
                this->update_Hessian_inverse(); 

                if(update_hessian_){
                    this->update_Hessian(); 
                }

                return true; 
            }

            virtual Matrix & get_Hessian() 
            {
                return H_prev_; 
            }

            virtual bool apply_Hinv(const Vector & g , Vector & s) const override
            {
                if(this->initialized())
                    s = H_prev_inv_ * g; 
                else
                    utopia_error("BFGS::apply_Hinv: Initialization needs to be done first. \n"); 
             
                return true; 
            }


            virtual Scalar compute_uHinvv_dot(const Vector & u, const Vector & v) const override
            {
                if(this->initialized())
                    return dot(u, H_prev_inv_ * v); 
                else{
                    utopia_error("BFGS::compute_uHinvv_dot: Initialization needs to be done first. \n"); 
                    return false; 
                }
            }


            virtual bool apply_H(const Vector & v, Vector & result) const override
            {
                if(update_hessian_ && this->initialized())
                    result = H_prev_ * v; 
                else
                   utopia_error("BFGS::apply_H can be used only, if H is computed. \n Please turn on update_hessian option. \n"); 

                return true; 
            }


            virtual Scalar compute_uHv_dot(const Vector & u, const Vector & v) const override
            {
                if(update_hessian_ && this->initialized())
                    return dot(u, H_prev_ * v); 
                else
                {
                    utopia_error("BFGS::compute_uHv_dot can be used only, if H is computed. \n Please turn on update_hessian option. \n"); 
                    return 0; 
                }
            }

            virtual Scalar compute_uHu_dot(const Vector & u) const override
            {
                if(update_hessian_ && this->initialized())
                    return dot(u, H_prev_ * u);
                else
                { 
                    utopia_error("BFGS::compute_uHu_dot can be used only, if H is computed. \n Please turn on update_hessian option. \n"); 
                    return 0; 
                }
            }

            void set_update_hessian(const bool flg)
            {
                update_hessian_ = flg; 
            }
            
            bool get_update_hessian() const
            {
                return update_hessian_; 
            }


        private:
            virtual void update_Hessian_inverse()
            {
                Scalar s_y =  dot(s_, y_); 
                Scalar yHy = dot(y_, H_prev_inv_ * y_);

                // checking for numerical instabilities
                if(s_y < this->num_tol() || yHy < this->num_tol() || !std::isfinite(s_y) || !std::isfinite(yHy))
                {
                    // std::cout<<"--- BFGS::update_Hessian_inverse is reaching num. tolerance \n"; 
                    return;
                }

                Matrix ss  = outer(s_, s_);
                Vector Hy  = H_prev_inv_ * y_; 
                Matrix Hys = outer(Hy, s_); 
                
                Matrix sy_outerH  = outer(s_, y_);
                sy_outerH = sy_outerH *  H_prev_inv_; 

                H_prev_inv_ += (s_y + yHy)/(s_y*s_y) * ss; 
                H_prev_inv_ -= (1./s_y) * (Hys + sy_outerH);  
            }


            virtual void update_Hessian()
            {
                Scalar y_s = dot(y_, s_);

                Vector H_s = H_prev_ * s_;
                Scalar sHs = dot(s_, H_s);
                
                if(y_s < this->num_tol() || !std::isfinite(y_s) || sHs < this->num_tol() || !std::isfinite(sHs)){
                    // std::cout<<"--- BFGS::update_Hessian is reaching num. tolerance \n"; 
                    return;
                }

                Matrix yy = outer(y_, y_); 
                H_prev_ += 1./y_s * yy; 

                Matrix a  = outer(H_s, H_s);
                H_prev_ -= 1./sHs * a;
            }


        private:
            Vector s_;              // x_{k+1} - x_{k}
            Vector y_;              // g_{k+1} - g_{k}
            Matrix H_prev_inv_;     // H^{-1}_{k}
            Matrix H_prev_;         // H^{1}_{k}
            bool update_hessian_; 

        };

}

#endif //UTOPIA_HESSIAN_BFGS_HPP
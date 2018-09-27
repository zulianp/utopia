#ifndef UTOPIA_HESSIAN_APPROX_HPP
#define UTOPIA_HESSIAN_APPROX_HPP

#include "utopia_Core.hpp"


namespace utopia
{

        template<class Matrix, class Vector>
        class HessianApproximation
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


            public:

                HessianApproximation(): num_tol_(1e-12), initialized_(false)
                {

                }

                virtual ~HessianApproximation() { }

                virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x) = 0; 

                Scalar num_tol(){return num_tol_; } 
                void num_tol(Scalar & tol ){num_tol_ = tol; }

                void initialized(const bool init){initialized_ = init; }
                bool initialized(){return initialized_; }

                // applications of inverse of Hessian
                virtual bool apply_Hinv(const Vector & /* g */, Vector & /*s */){return false; }
                virtual Scalar compute_uHinvv_dot(const Vector &/*u*/, const Vector & /*v*/){return false; }
                
                // applications of Hessian
                virtual bool apply_H(const Vector & /*v*/ , Vector & /*r */){return false; }
                virtual Scalar compute_uHv_dot(const Vector &/*u*/, const Vector & /*v*/){return false; }
                virtual Scalar compute_uHu_dot(const Vector &/*u*/){return false; }

                virtual bool update(const Vector & /* s  */, const Vector &  /* y */ ) = 0; 

                // to be deleted at some point
                virtual bool approximate_hessian( Function<Matrix, Vector> &/*fun*/,  const Vector & /*sol_new*/, const Vector & /*step*/, Matrix & /*hessian_old_new*/,  Vector & /*grad_old_new*/) {return 0; }
                virtual bool approximate_hessian_inverse( Function<Matrix, Vector> & /*fun*/,  const Vector & /*sol_new*/, const Vector & /*step*/, Matrix & /*hessian_old_new*/,  Vector & /*grad_old_new*/) {return 0; }


            private:
                Scalar num_tol_;
                bool initialized_; 

        };


        template<class Matrix, class Vector>
        class BFGS : public HessianApproximation<Matrix, Vector>
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

            public:

                BFGS(): HessianApproximation<Matrix, Vector>(), update_hessian_(false)
                {

                }

                virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x) override
                {
                    SizeType n = local_size(x).get(0);

                    if(update_hessian_)
                    {
                        if(!fun.hessian(x, H_prev_))
                            H_prev_ = local_identity(n, n);
                    }

                    H_prev_inv_ = local_identity(n, n);

                    s_ = local_zeros(n); 
                    y_ = local_zeros(n); 

                    this->initialized(true); 

                    return true;
                }                

                virtual bool update(const Vector & s_in, const Vector & y_in ) override
                {
                    s_ = s_in; 
                    y_ = y_in;
                    
                    this->update_Hessian_inverse(); 

                    if(update_hessian_){
                        this->update_Hessian(); 
                    }

                    return true; 
                }

                virtual bool apply_Hinv(const Vector & g , Vector & s) override
                {
                    if(this->initialized())
                        s = H_prev_inv_ * g; 
                    else
                        utopia_error("BFGS::apply_Hinv: Initialization needs to be done first. \n"); 
                 
                    return true; 
                }


                virtual Scalar compute_uHinvv_dot(const Vector & u, const Vector & v) override
                {
                    if(this->initialized())
                        return dot(u, H_prev_inv_ * v); 
                    else{
                        utopia_error("BFGS::compute_uHinvv_dot: Initialization needs to be done first. \n"); 
                        return false; 
                    }
                }


                virtual bool apply_H(const Vector & v, Vector & result) override
                {
                    std::cout<<"apply_H: update_hessian_: "<< update_hessian_ << "  \n"; 
                    if(update_hessian_ && this->initialized())
                        result = H_prev_ * v; 
                    else
                       utopia_error("BFGS::apply_H can be used only, if H is computed. \n Please turn on update_hessian option. \n"); 

                    return true; 
                }


                virtual Scalar compute_uHv_dot(const Vector & u, const Vector & v) override
                {
                    if(update_hessian_ && this->initialized())
                        return dot(u, H_prev_inv_ * v); 
                    else
                    {
                        utopia_error("BFGS::compute_uHv_dot can be used only, if H is computed. \n Please turn on update_hessian option. \n"); 
                        return 0; 
                    }
                }

                virtual Scalar compute_uHu_dot(const Vector & u) override
                {
                    if(update_hessian_ && this->initialized())
                        return dot(u, H_prev_inv_ * u);
                    else
                    { 
                        utopia_error("BFGS::compute_uHu_dot can be used only, if H is computed. \n Please turn on update_hessian option. \n"); 
                        return 0; 
                    }
                }

                void set_update_hessian(const bool flg){update_hessian_ = flg; }
                bool get_update_hessian(){return update_hessian_; }

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
                static_assert(!utopia::is_sparse<Matrix>::value, "BFGS does not support sparse matrices."); 
                Vector s_;              // x_{k+1} - x_{k}
                Vector y_;              // g_{k+1} - g_{k}
                Matrix H_prev_inv_;     // H^{-1}_{k}
                Matrix H_prev_;         // H^{1}_{k}
                bool update_hessian_; 

            };


        // template<class Matrix, class Vector>
        // class SR1 : public HessianApproximation<Matrix, Vector>
        // {
        //     typedef UTOPIA_SCALAR(Vector)    Scalar;
        //     typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        //     public:

        //         virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x) override
        //         {
        //             SizeType n = local_size(x).get(0);

        //             if(!fun.hessian(x, H_prev_))
        //                 H_prev_ = local_identity(n, n);

        //             s_ = local_zeros(n); 
        //             y_ = local_zeros(n); 

        //             return true;
        //         }                    

        //         virtual bool approximate_hessian(
        //             Function<Matrix, Vector> &fun,
        //             const Vector & x_new,
        //             const Vector & s,
        //             Matrix &hessian_old_new,
        //             Vector &grad_old_new) override
        //         {
        //             Vector g_diff = grad_old_new;
        //             if(!fun.gradient(x_new, grad_old_new)) return false;
        //             g_diff = grad_old_new - g_diff;

        //             Vector help_vec = g_diff - hessian_old_new * s;
        //             Scalar denom = dot(help_vec, s);

        //             // here, we prevent numerical instabilities leading to very small denominator
        //             if(std::abs(denom) >= this->num_tol() * norm2(s) * norm2(help_vec))
        //                     hessian_old_new += 1./denom * Matrix(outer(help_vec, help_vec));
        //             return true;
        //         }

        //         virtual bool update(const Vector & s_in, const Vector & y_in ) override
        //         {
        //             s_ = s_in; 
        //             y_ = y_in;
        //             return true; 
        //         }


        //     private:
        //         static_assert(!utopia::is_sparse<Matrix>::value, "SR1 does not support sparse matrices."); 
        //         Vector s_;  // x_{k+1} - x_{k}
        //         Vector y_;  // g_{k+1} - g_{k}
        //         Matrix H_prev_; 

        // };

}

#endif //UTOPIA_HESSIAN_APPROX_HPP
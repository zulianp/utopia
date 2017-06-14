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
            virtual ~HessianApproximation() { }

            virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x, Matrix &H)
            {
                // (void) fun;
                // SizeType n_local = local_size(x).get(0);
                // H = local_identity(n_local, n_local);

                //for sparse
                fun.hessian(x, H); 
                
                // opt 1
                // fun.hessian(x, H); 
                // H = 0.0 * H + Matrix(diag(diag(H)));

                //opt 2
                //hessian *= 0.;
                //hessian += local_identity(n_local, n_local);
                
                return true; 
            }



            virtual bool approximate_hessian( Function<Matrix, Vector> &fun,  const Vector &sol_new, const Vector &step, Matrix &hessian_old_new,  Vector &grad_old_new) = 0;
            virtual bool approximate_hessian_inverse( Function<Matrix, Vector> &fun,  const Vector &sol_new, const Vector &step, Matrix &hessian_old_new,  Vector &grad_old_new) {return 0; };
        };
        

        template<class Matrix, class Vector>
        class BFGS : public HessianApproximation<Matrix, Vector>
        {
            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

            public:
                virtual bool approximate_hessian(
                    Function<Matrix, Vector> &fun, 
                    const Vector &sol_new,
                    const Vector &step,
                    Matrix &hessian_old_new, 
                    Vector &grad_old_new) override
                {
                    Vector diff_grad = grad_old_new;
                    if(!fun.gradient(sol_new, grad_old_new)) return false;
                    diff_grad = grad_old_new - diff_grad;

                    Vector H_x_step = hessian_old_new * step;
                    hessian_old_new += (1./Scalar(dot(diff_grad, step)) * Matrix(outer(diff_grad, diff_grad)));

                    //utopia::is_sparse<Matrix>::value
                    // temp += (outer(H_x_step, H_x_step)/dot(step, H_x_step));
                    
                    Matrix a  = outer(H_x_step, H_x_step);
                    Scalar b = dot(step, H_x_step); 
                    hessian_old_new -= (1./b * a); 


                    return true;
                }   






            };



        template<class Matrix, class Vector>
        class SR1 : public HessianApproximation<Matrix, Vector>
        {

            typedef UTOPIA_SCALAR(Vector)    Scalar;
            typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

            public:

                SR1(): num_tol_(1e-8)
                { }

                virtual bool approximate_hessian(
                    Function<Matrix, Vector> &fun, 
                    const Vector & x_new,
                    const Vector & s,
                    Matrix &hessian_old_new, 
                    Vector &grad_old_new) override
                {
                    Vector g_diff = grad_old_new;
                    if(!fun.gradient(x_new, grad_old_new)) return false;
                    g_diff = grad_old_new - g_diff;

                    Vector help_vec = g_diff - hessian_old_new * s; 
                    Scalar denom = dot(help_vec, s); 

                    // here, we prevent numerical instabilities leading to very small denominator 
                    if(std::abs(denom) >= num_tol_ * norm2(s) * norm2(help_vec))
                            hessian_old_new += 1./denom * Matrix(outer(help_vec, help_vec)); 
                    return true;
                }   


            private:
                Scalar num_tol_; 

            };

}

#endif //UTOPIA_HESSIAN_APPROX_HPP        
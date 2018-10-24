#ifndef UTOPIA_LBFGS_HPP
#define UTOPIA_LBFGS_HPP

#include "utopia_Core.hpp"


namespace utopia
{
    // TO BE DONE.... 
    template<class Matrix, class Vector>
    class LBFGS : public HessianApproximation<Matrix, Vector>
    {

        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        public:

            LBFGS(const SizeType & m): m_(m), current_m_(0), theta_(1.0), gamma_(1.0)
            {

            }

            virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x) override
            {
                theta_ = 1.0; 
                gamma_ = 1.0; 

                SizeType n = local_size(x).get(0);
                H0_ = local_identity(n, n);

                current_m_ = 0; 
                this->initialized(true); 

                Y_.resize(m_); 
                S_.resize(m_); 

                a_.resize(m_); 
                b_.resize(m_);                 

                rho_.resize(m_); 

                return true;
            }   

            virtual bool update(const Vector &  s  , const Vector &  y ) override
            {
                
                if(!this->initialized())
                {
                    utopia_error("BFGS::update: Initialization needs to be done before updating. \n"); 
                    return false; 
                }

                Scalar nom      = dot(y,y);
                Scalar denom    = dot(y,s); 

                // theta and gamma are inverse of each other
                theta_ = nom/denom; 
                gamma_ = denom/nom; 


                // if denom > eps, hessian approx. should be positive semidefinite
                if(denom < 1e-12 || !std::isfinite(denom))
                {
                    // if(mpi_world_rank()==0)
                    //     utopia_warning("L-BFGS-B: Curvature condition not satified. Skipping update. \n"); 

                    return false; 
                }

                if(m_ ==0)
                    return true; 

                if(current_m_ < m_)
                {
                    Y_[current_m_]      = y; 
                    S_[current_m_]      = s; 
                    rho_[current_m_]    = 1./denom; 

                    Scalar ys = 1./std::pow(denom, 1./2.); 
                    b_[current_m_] = ys * y; 

                }
                else
                {
                    Y_[0]   = y; 
                    S_[0]   = s; 
                    rho_[0] = 1./denom; 

                    Scalar ys = 1./std::pow(denom, 1./2.); 
                    b_[0] = ys * y;                     

                    std::rotate(Y_.begin(), Y_.begin() + 1, Y_.end());
                    std::rotate(S_.begin(), S_.begin() + 1, S_.end());  

                    std::rotate(rho_.begin(), rho_.begin() + 1, rho_.end());  
                    std::rotate(b_.begin(), b_.begin() + 1, b_.end());  
                }

                current_m_++; 

                // this is here just to save computation time Bv product 
                this->update_a_b(); 

                return true; 
            }

            virtual bool apply_Hinv(const Vector & g, Vector & q) const override
            {
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 
                std::vector<Scalar> alpha_inv(current_memory_size); 
                q = g; 

                for(auto i=current_memory_size-1; i >=0; i--)
                {
                    alpha_inv[i] = rho_[i] * dot(S_[i], q); 
                    q -=  alpha_inv[i] * Y_[i]; 
                }

                q = gamma_ * (H0_ * q); 


                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar betta_inv = rho_[i] * dot(Y_[i], q); 
                    q += (alpha_inv[i] - betta_inv)  * S_[i];  
                }

                return true; 
            }

            virtual bool apply_H(const Vector & v , Vector & result) const  override
            {
                result = (theta_ * H0_) * v; 
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 

                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar bv  = dot(b_[i], v); 
                    Scalar av  = dot(a_[i], v); 

                    result += (bv * b_[i]) - (av * a_[i]); 
                }
                return true; 
            }


            virtual Matrix & get_Hessian() override
            {
                std::cerr<<"--- not implemented yet---- \n"; 
                return H0_;
            } 

            void set_memory_size(const SizeType & m)
            {
                m_ = m; 
            }

            SizeType get_memory_size() const 
            {
                return m_; 
            }


            void update_a_b()
            {
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 

                for(auto k =0; k < current_memory_size; k++)
                {
                    a_[k] = (theta_*H0_) * S_[k]; 

                    for(auto i = 0; i < k; i++)
                    {
                        Scalar bS = dot(b_[i], S_[k]) ; 
                        Scalar aS = dot(a_[i], S_[k]); 

                        a_[k] += bS * b_[i]; 
                        a_[k] -= aS * a_[i]; 
                    }

                    Scalar scaling_factor = (std::pow(dot(S_[k], a_[k]), 1./2.)); 
                    a_[k] = 1.0/scaling_factor * a_[k]; 
                }
            }


        private:
            // static_assert(utopia::is_sparse<Matrix>::value, "BFGS does not support sparse matrices."); 
            // static_assert(!utopia::is_sparse<DenseMatrix>::value, "BFGS does not support sparse matrices."); 

            SizeType m_; // memory size 
            SizeType current_m_; // current amount of vectors in the memory 
            Scalar theta_; 
            Scalar gamma_; 
            
            Matrix H0_; 

            std::vector<Vector> Y_; 
            std::vector<Vector> S_; 


            std::vector<Vector > b_; 
            std::vector<Vector > a_; 


            std::vector<Scalar > rho_; 


        };



}

#endif //UTOPIA_LBFGS_HPP

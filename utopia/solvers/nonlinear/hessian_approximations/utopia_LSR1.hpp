#ifndef UTOPIA_LSR1_HPP
#define UTOPIA_LSR1_HPP

#include "utopia_Core.hpp"


namespace utopia
{
    template<class Vector>
    class LSR1 final: public HessianApproximation<Vector>
    {

        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        public:

            LSR1(const SizeType & m): m_(m), current_m_(0), theta_(1.0), gamma_(1.0)
            {

            }

            virtual void initialize(const SizeType & n) override
            {
                theta_ = 1.0; 
                gamma_ = 1.0; 

                current_m_ = 0; 

                Y_.resize(m_); 
                S_.resize(m_); 

                p_.resize(m_);
                p_inv_.resize(m_);

                this->initialized(true); 
            }   

            virtual bool update(const Vector &  s, const Vector &  y ) override
            {
                
                if(!this->initialized())
                {
                    utopia_error("BFGS::update: Initialization needs to be done before updating. \n"); 
                    return false; 
                }
                else if(m_ == 0)
                    return true; 
            
                // if denom > eps, hessian approx. should be positive semidefinite
                // if(denom < 1e-12 || !std::isfinite(denom))
                // {
                //     // if(mpi_world_rank()==0)
                //     //     utopia_warning("L-BFGS-B: Curvature condition not satified. Skipping update. \n"); 

                //     return false; 
                // }


                if(current_m_ < m_)
                {
                    Y_[current_m_]      = y; 
                    S_[current_m_]      = s; 
                }
                else
                {
                    Y_[0]   = y; 
                    S_[0]   = s; 

                    std::rotate(Y_.begin(), Y_.begin() + 1, Y_.end());
                    std::rotate(S_.begin(), S_.begin() + 1, S_.end());  
                }

                current_m_++; 

                precompute_p(); 
                precompute_p_inv(); 

                return true; 
            }

            virtual bool apply_Hinv(const Vector & v, Vector & result) const override
            {
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 
                result = gamma_ * v;

                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar scaling_factor = dot(p_inv_[i], v)/dot(p_inv_[i], Y_[i]); 
                    result += scaling_factor * p_inv_[i];
                }

                return true; 
            }

            virtual bool apply_H(const Vector & v, Vector & result) const  override
            {
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 
                result = theta_ * v;

                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar scaling_factor = dot(p_[i], v)/dot(p_[i], S_[i]); 
                    result += scaling_factor * p_[i];
                }

                return true; 
            }

            void set_memory_size(const SizeType & m)
            {
                m_ = m; 
            }

            SizeType get_memory_size() const 
            {
                return m_; 
            }


        private:
            void precompute_p()
            {
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 

                for(auto i=0; i < current_memory_size; i++)
                {
                    p_[i] = Y_[i] - theta_ * S_[i];
                    for(auto j=0; j < i; j++)
                    {
                        Scalar scaling_factor = dot(p_[j], S_[i])/dot(p_[j], S_[j]); 
                        p_[i] -= scaling_factor* p_[j];  
                    }
                }
            }


            void precompute_p_inv()
            {
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 

                for(auto i=0; i < current_memory_size; i++)
                {
                    p_inv_[i] = S_[i] - gamma_ * Y_[i];
                    for(auto j=0; j < i; j++)
                    {
                        Scalar scaling_factor = dot(p_inv_[j], Y_[i])/dot(p_inv_[j], Y_[j]); 
                        p_inv_[i] -= scaling_factor* p_inv_[j];  
                    }
                }
            }


        private:
            SizeType m_; // memory size 
            SizeType current_m_; // current amount of vectors in the memory 
            
            Scalar theta_; 
            Scalar gamma_; 

            std::vector<Vector> Y_; 
            std::vector<Vector> S_; 

            std::vector<Vector> p_; 
            std::vector<Vector> p_inv_; 

        };

}

#endif //UTOPIA_LBFGS_HPP

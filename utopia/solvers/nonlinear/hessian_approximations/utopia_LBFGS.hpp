#ifndef UTOPIA_LBFGS_HPP
#define UTOPIA_LBFGS_HPP

#include "utopia_Core.hpp"


namespace utopia
{
    // TODO:: speedup by preallocating just for H_inv, or H products... 
    template<class Vector>
    class LBFGS final: public HessianApproximation<Vector>
    {

        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        public:

            LBFGS(const SizeType & m): m_(m), current_m_(0), theta_(1.0), gamma_(1.0)
            {
               
            }

            virtual void initialize() override
            {
                theta_ = 1.0; 
                gamma_ = 1.0; 

                current_m_ = 0; 

                Y_.resize(m_); 
                S_.resize(m_); 

                a_.resize(m_); 
                b_.resize(m_);                 

                rho_.resize(m_); 
                Sb_dots_.resize(m_, std::vector<Scalar>(m_));

                this->initialized(true); 
            }   


            virtual void reset() override
            {
                Y_.clear(); 
                S_.clear(); 

                a_.clear(); 
                b_.clear(); 

                rho_.clear(); 
                Sb_dots_.clear(); 

                this->initialize(); 
            }

            inline LBFGS<Vector> * clone() const override
            {
                return new LBFGS<Vector>(*this);
            }

            virtual bool update(const Vector &  s, const Vector &  y ) override
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

                    if(current_m_ > 0)
                        dots(S_[current_m_], b_, Sb_dots_[current_m_]); 
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

                    dots(S_[m_-1], b_, Sb_dots_[0]); 
                    std::rotate(Sb_dots_.begin(), Sb_dots_.begin() + 1, Sb_dots_.end());  

                    for(auto i=0; i < Sb_dots_.size()-1; i++){
                        for(auto j=1; j < Sb_dots_.size(); j++){
                            Sb_dots_[i][j-1] = Sb_dots_[i][j]; 
                        }
                    }

                }

                current_m_++; 

                // this is here just to save computation time Bv product 
                this->update_a_b(); 

                return true; 
            }

            virtual bool apply_Hinv(const Vector & g, Vector & q) const override
            {
                if(!this->initialized())
                    utopia_error("utopia::LBFGS::apply_Hinv:: missing initialization... \n"); 

                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 
                std::vector<Scalar> alpha_inv(current_memory_size); 
                q = g; 

                for(auto i=current_memory_size-1; i >=0; i--)
                {
                    alpha_inv[i] = rho_[i] * dot(S_[i], q); 
                    q -=  alpha_inv[i] * Y_[i]; 
                }

                q = gamma_ * q; 

                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar betta_inv = rho_[i] * dot(Y_[i], q); 
                    q += (alpha_inv[i] - betta_inv)  * S_[i];  
                }

                if(has_nan_or_inf(q))
                    q = gamma_ * g;

                return true; 
            }

            virtual bool apply_H(const Vector & v , Vector & result) const  override
            {
                if(!this->initialized())
                    utopia_error("utopia::LBFGS::apply_Hinv:: missing initialization... \n"); 
                
                result = theta_ * v; 
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 

                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar bv  = dot(b_[i], v); 
                    Scalar av  = dot(a_[i], v); 

                    result += (bv * b_[i]) - (av * a_[i]); 
                }

                if(has_nan_or_inf(result))
                    result = theta_ * v;

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
            void update_a_b()
            {
                if(!this->initialized())
                    utopia_error("utopia::LBFGS::apply_Hinv:: missing initialization... \n"); 

                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_; 

                for(auto k =0; k < current_memory_size; k++)
                {
                    a_[k] = theta_ * S_[k]; 

                    for(auto i = 0; i < k; i++)
                    {
                        Scalar aS = dot(a_[i], S_[k]); 

                        a_[k] += Sb_dots_[k][i] * b_[i]; 
                        a_[k] -= aS * a_[i]; 
                    }

                    Scalar scaling_factor = (std::pow(dot(S_[k], a_[k]), 1./2.)); 
                    a_[k] = 1.0/scaling_factor * a_[k]; 
                }
            }


        private:
            SizeType m_; // memory size 
            SizeType current_m_; // current amount of vectors in the memory 
            Scalar theta_; 
            Scalar gamma_; 
            
            std::vector<Vector> Y_; 
            std::vector<Vector> S_; 

            std::vector<Vector > b_; 
            std::vector<Vector > a_; 
            std::vector<Scalar > rho_; 

            std::vector<std::vector<Scalar> > Sb_dots_; 


        };

}

#endif //UTOPIA_LBFGS_HPP

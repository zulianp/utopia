#ifndef UTOPIA_LBFGS_HPP
#define UTOPIA_LBFGS_HPP

#include "utopia_Core.hpp"
#include "utopia_HessianApproximation.hpp"


namespace utopia
{
    enum LBFGSDampingTechnique  {   NOCEDAL = 0,  
                                    POWEL   = 1};

    enum LBFGSScalingTechnique  {   NONE        = 0,  
                                    INITIAL     = 1,
                                    ADAPTIVE    = 2,
                                    FORBENIUS   = 3 };                                


    // TODO:: speedup by preallocating just for H_inv, or H products...
    template<class Vector>
    class LBFGS final: public HessianApproximation<Vector>
    {

        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;


        public:

            LBFGS(const SizeType & m):  m_(m), 
                                        current_m_(0), 
                                        theta_(1.0), 
                                        gamma_(1.0), 
                                        theta_min_(1.0), 
                                        skip_update_treshold_(1e-12), 
                                        damping_tech_(NOCEDAL), 
                                        scaling_tech_(ADAPTIVE)
            {

            }

            void initialize(const Vector & /*x_k*/, const Vector & /* g */) override
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


            void reset() override
            {
                Y_.clear();
                S_.clear();

                a_.clear();
                b_.clear();

                rho_.clear();
                Sb_dots_.clear();

                Vector x, g; 
                this->initialize(x, g);
            }

            inline LBFGS<Vector> * clone() const override
            {
                return new LBFGS<Vector>(*this);
            }

            bool update(const Vector &  s, const Vector &  y, const Vector & /*x*/, const Vector &  /* g */ ) override
            {

                if(!this->initialized())
                {
                    utopia_error("BFGS::update: Initialization needs to be done before updating. \n");
                    return false;
                }

                if(m_ ==0){
                    return true;
                }


                Vector y_hat = y; 
                bool skip_update = init_damping(y, s, y_hat);
                
                if(skip_update)
                {
                    this->init_scaling_factors(y_hat, s); 
                    return true; 
                }

                Scalar denom    = dot(y_hat,s);

                if(current_m_ < m_)
                {
                    Y_[current_m_]      = y_hat;
                    S_[current_m_]      = s;
                    rho_[current_m_]    = 1./denom;

                    Scalar ys = 1./std::pow(denom, 1./2.);
                    b_[current_m_] = ys * y_hat;

                    if(current_m_ > 0){
                        dots(S_[current_m_], b_, Sb_dots_[current_m_]);
                    }
                }
                else
                {
                    Y_[0]   = y_hat;
                    S_[0]   = s;
                    rho_[0] = 1./denom;

                    Scalar ys = 1./std::pow(denom, 1./2.);
                    b_[0] = ys * y_hat;

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

                // updating the factors... 
                this->init_scaling_factors(y_hat, s); 

                return true;
            }

            bool apply_Hinv(const Vector & g, Vector & q) const override
            {
                if(!this->initialized()){
                    utopia_error("utopia::LBFGS::apply_Hinv:: missing initialization... \n");
                }

                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;
                std::vector<Scalar> alpha_inv(current_memory_size);
                q = g;

                for(auto i=current_memory_size-1; i >=0; i--)
                {
                    alpha_inv[i] = rho_[i] * dot(S_[i], q);
                    q -=  alpha_inv[i] * Y_[i];
                }

                Vector q2 = q;
                this->apply_H0_inv(q2, q); 

                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar betta_inv = rho_[i] * dot(Y_[i], q);
                    q += (alpha_inv[i] - betta_inv)  * S_[i];
                }

                if(has_nan_or_inf(q)){
                    this->apply_H0_inv(g, q); 
                }

                return true;
            }

            bool apply_H(const Vector & v , Vector & result) const  override
            {
                // this->apply_H_fast(v, result); 

                // this->apply_H_slow(v, result); 

                this-apply_H_new(v, result); 

                return true;
            }


            bool apply_H_fast(const Vector & v , Vector & result) const
            {
                if(!this->initialized()){
                    utopia_error("utopia::LBFGS::apply_H:: missing initialization... \n");
                }

                result = theta_ * v;
                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar bv  = dot(b_[i], v);
                    Scalar av  = dot(a_[i], v);

                    result += (bv * b_[i]) - (av * a_[i]);
                }

                if(has_nan_or_inf(result)){
                    result = theta_ * v;
                }

                return true; 
            }



            bool apply_H_slow(const Vector & v , Vector & result) const
            {
                if(!this->initialized()){
                    utopia_error("utopia::LBFGS::apply_H:: missing initialization... \n");
                }


                std::vector<Vector> a, b; 
                a.resize(m_);
                b.resize(m_);

                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

                for(auto k =0; k < current_memory_size; k++)
                {
                    Scalar mult_b = 1./std::sqrt(dot(Y_[k], S_[k]));
                    b[k] = mult_b * Y_[k];

                    // a[k] = theta_ * S_[k]; 
                    this->apply_H0(S_[k], a[k]); 

                    for(auto i = 0; i < k; i++)
                    {
                        Scalar bS = dot(b[i], S_[k]);
                        Scalar aS = dot(a[i], S_[k]);

                        a[k] += bS * b[i];
                        a[k] -= aS * a[i];
                    }

                    Scalar mult_a = 1./std::sqrt(dot(S_[k], a[k]));
                    a[k] = mult_a * a[k]; 
                }                

                // result = theta_ * v;
                this->apply_H0(v, result); 

                for(auto i=0; i < current_memory_size; i++)
                {
                    Scalar bv  = dot(b[i], v);
                    Scalar av  = dot(a[i], v);

                    result += (bv * b[i]);
                    result -= (av * a[i]);
                }


                return true;
            }

            bool apply_H_new(const Vector & v , Vector & result) const
            {
                if(!this->initialized()){
                    utopia_error("utopia::LBFGS::apply_H:: missing initialization... \n");
                }

                std::vector<Vector> P; 
                std::vector<Scalar> stp, yts;
                stp.resize(m_);
                yts.resize(m_);
                P.resize(m_);

                SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

                for(auto i =0; i < current_memory_size; i++)
                {
                    this->apply_H0(S_[i], P[i]); 

                    yts[i] = dot(Y_[i], S_[i]); 

                    for(auto j=0; j < i; j++)
                    {
                        Scalar mult1 = dot(S_[j], P[i])/stp[j]; 
                        Scalar mult2 = dot(S_[i], Y_[j])/yts[j]; 

                        P[i] = P[i]  - (mult1 * P[j]) + (mult2 * Y_[j]); 
                    }

                    stp[i] = dot(S_[i], P[i]); 
                }                

                this->apply_H0(v, result); 

                for(auto i =0; i < current_memory_size; i++)
                {
                    Scalar stz = dot(S_[i], result); 
                    Scalar ytx = dot(Y_[i], v); 

                    result = result + ((ytx/yts[i]) * Y_[i]) - ((stz/stp[i]) * P[i]);
                }

                return true;
            }




            void memory_size(const SizeType & m)
            {
                m_ = m;
            }

            SizeType memory_size() const
            {
                return m_;
            }


            Scalar theta_min() const
            {
                return theta_min_; 
            }


            LBFGSDampingTechnique damping_tech() const
            {
                return damping_tech_; 
            }

            void damping_tech(const LBFGSDampingTechnique & damping)
            {
                damping_tech_ = damping; 
            }


            LBFGSScalingTechnique scaling_tech() const
            {
                return scaling_tech_; 
            }

            void scaling_tech(const LBFGSScalingTechnique & scaling)
            {
                scaling_tech_ = scaling; 
            }


            void theta_min(const Scalar & theta_min)
            {
                theta_min_ = theta_min; 
            }

            void read(Input &in) override
            {
                HessianApproximation<Vector>::read(in);
                in.get("memory_size", m_);
            }

            void print_usage(std::ostream &os) const override
            {
                HessianApproximation<Vector>::print_usage(os);
                this->print_param_usage(os, "memory_size", "int", "History/Memory size.", "5");
            }


            void H0_action(const std::function< void(const Vector &, Vector &) > operator_action)
            {
                H0_action_ = operator_action; 
            }

            void H0_inv_action(const std::function< void(const Vector &, Vector &) > operator_action)
            {
                H0_inv_action_ = operator_action; 
            }            


        private:
            void update_a_b()
            {
                if(!this->initialized()){
                    utopia_error("utopia::LBFGS::apply_Hinv:: missing initialization... \n");
                }

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



            bool init_damping(const Vector & y, const Vector & s, Vector & y_hat)
            {
                if(damping_tech_ == LBFGSDampingTechnique::POWEL)
                {
                    return init_damping_powel(y, s, y_hat); 
                }
                else
                {
                    return init_damping_nocedal(y, s, y_hat); 
                }
            }


            void init_scaling_factors(const Vector & y, const Vector &s)
            {
                if(scaling_tech_ == LBFGSScalingTechnique::NONE)
                {
                    theta_ = 1.0; 
                    gamma_ = 1.0; 
                    D_ = local_values(local_size(y).get(0), 1.0); 
                    D_inv_ = local_values(local_size(y).get(0), 1.0); 
                }
                else if(scaling_tech_ == LBFGSScalingTechnique::INITIAL)
                {
                    if(current_m_ ==0 )
                    {   
                        theta_ = dot(y, y)/ dot(y,s);
                        theta_ = std::max(theta_min_, theta_);
                        gamma_ = 1./theta_;           

                        if(!std::isfinite(theta_) || !std::isfinite(gamma_))
                        {
                            theta_ = 1.0; 
                            gamma_ = 1.0; 
                        }                     

                        D_ = local_values(local_size(y).get(0), theta_); 
                        D_inv_ = local_values(local_size(y).get(0), gamma_); 
                    }
                }
                else if(scaling_tech_ == LBFGSScalingTechnique::ADAPTIVE)
                {
                    theta_ = dot(y, y)/ dot(y,s);
                    theta_ = std::max(theta_min_, theta_);
                    gamma_ = 1./theta_;           

                    if(!std::isfinite(theta_) || !std::isfinite(gamma_))
                    {
                        theta_ = 1.0; 
                        gamma_ = 1.0; 
                    }                     

                    D_ = local_values(local_size(y).get(0), theta_); 
                    D_inv_ = local_values(local_size(y).get(0), gamma_); 

                }
                // FORBENIUS
                else
                {
                    SizeType current_memory_size = (current_m_ < m_) ? current_m_ : m_;

                    // TBD:: reallocate
                    std::vector<Vector> SY_point_wise_; 
                    std::vector<Vector> YY_point_wise_; 

                    SY_point_wise_.resize(m_); 
                    YY_point_wise_.resize(m_); 

                    for(auto i=0; i <current_memory_size; i++)
                    {
                        SY_point_wise_[i] = e_mul(S_[i], Y_[i]); 
                        YY_point_wise_[i] = e_mul(Y_[i], Y_[i]); 
                    }

                    Vector SY_sum = local_values(local_size(y).get(0), 0.0); 
                    Vector YY_sum = local_values(local_size(y).get(0), 0.0); 
                

                    for(auto i=0; i <current_memory_size; i++)
                    {
                        SY_sum += SY_point_wise_[i];  
                        YY_sum += YY_point_wise_[i]; 
                    }

                    Scalar theta_working = dot(y, y)/ dot(y,s);
                    theta_working = std::max(theta_min_, theta_working);       

                    if(!std::isfinite(theta_working))
                    {
                        theta_working = 1.0; 
                    }     

                    Vector D_help =  0.0*y; 
                    Vector D_inv_help = 0.0*y; 

 
                    {
                        Read<Vector> r_ub(YY_sum), r_lb(SY_sum);
                        Write<Vector> wv(D_help);

                        each_write(D_help, [&YY_sum, &SY_sum, theta_working](const SizeType i) -> double 
                        {
                            Scalar yy =  YY_sum.get(i); Scalar sy =  SY_sum.get(i); 

                            if(!std::isfinite(yy) || !std::isfinite(sy))
                            {
                                return theta_working; 
                            }

                            Scalar dii = (yy/sy); 

                            if(dii > 10e-10 && std::isfinite(dii))
                            {
                                return dii; 
                            }
                            else
                            {
                                return theta_working;
                            }

                        }   );   
                    }

                    {
                        Read<Vector> r_ub(D_help); 
                        Write<Vector> wv(D_inv_help);

                        each_write(D_inv_help, [&D_help, theta_working](const SizeType i) -> double 
                        {
                            Scalar dii =  D_help.get(i);

                            if(std::isfinite(dii) && dii!=0.0)
                            {
                                return 1./dii; 
                            }
                            else
                            {
                                return theta_working; 
                            }

                        } );   
                    }

                    D_ = D_help; 
                    D_inv_ = D_inv_help; 
                }

            }


            bool init_damping_nocedal(const Vector & y, const Vector & s, Vector & y_hat)
            {
                bool skip_update = false; 
                if(dot(y,s) < skip_update_treshold_)
                {
                    if(mpi_world_rank()==0){
                        utopia_warning("L-BFGS-B: Curvature condition not satified. Skipping update. \n");
                    }

                    skip_update = true; 
                }                

                y_hat = y; 

                return skip_update; 
            }



            bool init_damping_powel(const Vector & y, const Vector & s, Vector & y_hat)
            {
                Vector Bs = 0.0*s; 
                this->apply_H(s, Bs); 

                Scalar sy   = dot(y,s);
                Scalar sBs  = dot(s, Bs); 

                Scalar psi; 

                if(sy >= 0.2 * sBs){
                    psi = 1.0;
                }
                else{
                    psi = (0.8* sBs)/(sBs - sy); 
                }

                y_hat = (psi*y) +((1.0-psi) * Bs);

                return false; 
            }


            void apply_H0(const Vector & v, Vector & result) const 
            {

                if(H0_action_)
                {
                    H0_action_(v, result); 
                }
                else
                {
                    if(current_m_ ==0)
                    {
                        result = theta_ * v; 
                    }
                    else
                    {
                        result = e_mul(D_, v); 
                    }
                }
            }

            void apply_H0_inv(const Vector & v, Vector & result) const 
            {

                if(H0_inv_action_)
                {
                    H0_inv_action_(v, result); 
                }
                else
                {
                    if(current_m_ ==0)
                    {
                        result = gamma_ * v; 
                    }
                    else
                    {
                        result = e_mul(D_inv_, v); 
                    }
                }
            }            


        private:
            SizeType m_; // memory size
            SizeType current_m_; // current amount of vectors in the memory

            Scalar theta_;
            Scalar gamma_;
            Scalar theta_min_; 
            Scalar skip_update_treshold_; 

            std::vector<Vector> Y_;
            std::vector<Vector> S_;

            std::vector<Vector > b_;
            std::vector<Vector > a_;
            std::vector<Scalar > rho_;

            std::vector<std::vector<Scalar> > Sb_dots_;

            LBFGSDampingTechnique damping_tech_; 
            LBFGSScalingTechnique scaling_tech_; 

            Vector D_; 
            Vector D_inv_; 

            // std::vector<Vector> SY_point_wise_; 
            // std::vector<Vector> YY_point_wise_; 


            std::function< void(const Vector &, Vector &) > H0_action_; 
            std::function< void(const Vector &, Vector &) > H0_inv_action_;            

        };


// TBD:: - optimize dot products, fix pre-allocations, fix dots 

}

#endif //UTOPIA_LBFGS_HPP

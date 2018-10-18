#ifndef UTOPIA_LBFGSB_HPP
#define UTOPIA_LBFGSB_HPP

#include "utopia_Core.hpp"

namespace utopia
{

template<class Matrix, class DenseMatrix, class Vector>
class LBFGSB : public HessianApproximation<Matrix, Vector>
{

    typedef UTOPIA_SCALAR(Vector)    Scalar;
    typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

    typedef utopia::LinearSolver<Matrix, Vector>    LinSolver;

    public:

        LBFGSB( const SizeType & m, 
                const std::shared_ptr <LinSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >()):
                m_(m), current_m_(0), linear_solver_(linear_solver)
        {
            // to be factored out.... 
            cp_.set_memory_size(-1); 
        }


        virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x) override
        {
            SizeType n = local_size(x).get(0);
            H_ = local_identity(n, n);


            if(cp_.get_memory_size()==-1)
                cp_.set_memory_size(size(x).get(0));

            this->initialized(true);
            current_m_ = 0; 

            SizeType local_col_W = 0; 
            if(mpi_world_rank()==0)
                local_col_W = 2.0 * m_; 

            // W \in R^{n \times 2m }
            W_ = local_values(n, local_col_W, 0.0); 

            // M \in R^{2m \times 2m }
            M_ = local_values(local_col_W, local_col_W, 0.0); 

            SizeType local_col_YS = 0; 
            if(mpi_world_rank()==0)
                local_col_YS = m_; 


            Y_ = local_values(n, local_col_YS, 0.0); 
            S_ = local_values(n, local_col_YS, 0.0); 

            theta_ = 1.0; 

            d_elements_.resize(m_); 
            L_dots_.resize(m_, std::vector<Scalar>(m_));

            return true;
        }


        virtual void set_linear_solver(const std::shared_ptr<LinSolver> &linear_solver)
        {
            linear_solver_ = linear_solver; 
        }


        virtual bool update(const Vector &  s, const Vector &  y) override
        {

            if(!this->initialized())
            {
                utopia_error("BFGS::update: Initialization needs to be done before updating. \n"); 
                return false; 
            }

            Scalar denom    = dot(y,s); 
            Scalar nom      = dot(y,y);

            // if denom > eps, hessian approx. should be positive semidefinite
            if(denom < 1e-12)
            {
                // if(mpi_world_rank()==0)
                //     utopia_warning("L-BFGS-B: Curvature condition not satified. Skipping update. \n"); 

                return false; 
            }

            theta_ = nom/denom; 

            this->shift_cols_left_replace_last_col(Y_, y); 
            this->shift_cols_left_replace_last_col(S_, s); 

            this->buildL_d_elems(s); 

            this->buildW(); 
            this->buildM(); 

            current_m_++; 

            return true;
        }

        virtual bool apply_Hinv(const Vector & g, Vector & s) const override
        {
            this->apply_inverse_to_vec(g, W_, s); 
            return true;
        }


        virtual bool apply_reduced_Hinv(const Vector & feasible_set, const Vector & g, Vector & s) const
        {
            SizeType local_feasible_size = reduced_primal_method_.get_local_size_feasible_set(feasible_set); 
            Vector g_reduced = local_zeros(local_feasible_size);         
            reduced_primal_method_.build_reduced_vector(g, feasible_set, g_reduced); 

            DenseMatrix W_reduced = local_values(local_feasible_size, local_size(W_).get(1), 0.0);
            reduced_primal_method_.build_reduced_matrix(W_, feasible_set, W_reduced); 

            this->apply_inverse_to_vec(g_reduced, W_reduced, s); 

            return true;
        }


        virtual bool apply_H(const Vector & v, Vector & result) const override
        {
            Vector Y_v = transpose(Y_) * v;
            Vector S_v = theta_ * transpose(S_) * v;

            Vector p = Vector(Blocks<Vector>(
            {
                make_ref(Y_v), make_ref(S_v)
            }));

            // if(current_m_ <= m_)
            // {
                this->apply_M(p, result); 
            // }
            // else
            // {
            //     result = M_*p; 
            // }


            result = (theta_ * v) - W_*result; 

            return false; 
        }


        virtual Matrix & get_Hessian() override
        {
            if(current_m_ > 1)
            {
                H_ = local_identity(local_size(H_).get(0), local_size(H_).get(1)); 
                H_ = theta_ * H_; 

                DenseMatrix Y_T = transpose(Y_); 
                DenseMatrix S_T = theta_ * transpose(S_); 

                DenseMatrix P =  DenseMatrix(Blocks<DenseMatrix>(2, 1,
                {
                    make_ref(Y_T), make_ref(S_T)
                }));

                DenseMatrix result; 
                // if(current_m_ <= m_)
                // {
                    this->apply_M(P, result);    
                // }
                // else
                // {
                //     result = M_*P; 
                // }

                H_ =  H_ - (transpose(P) * result);     
            }
            else
            {
                H_ = local_identity(local_size(H_).get(0), local_size(H_).get(1)); 
            }
            
            // utopia_warning("LBFGS::get_Hessian returns dense matrix ...."); 

            return H_; 
        }

        virtual Matrix & get_Hessian_inv() override
        {
            std::cerr << "--- not implemented yet---- \n";
            return H_;
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
        void shift_cols_left_replace_last_col(DenseMatrix & M, const Vector & col) const
        {
            const DenseMatrix M_old = M;

            {  
                Write<DenseMatrix>  write(M);

                Read<DenseMatrix>   read(M_old);
                Read<Vector>      read_vec(col);

                Range r = row_range(M_old);
                Range c = col_range(M_old);

                for (SizeType i = r.begin(); i != r.end(); ++i) {
                    for (SizeType j = c.begin(); j != c.end() - 1; ++j) {
                        M.set(i, j, M_old.get(i, j + 1));
                    }
                    M.set(i, c.end() - 1, col.get(i));
                }
            } 
        }

        void buildW()
        {
            DenseMatrix S_theta = theta_ * S_; 
            W_ = DenseMatrix(Blocks<DenseMatrix>(1, 2,
            {
                make_ref(Y_), make_ref(S_theta)
            }));
        }


        void buildL_d_elems(const Vector &s)
        {
            Vector dots = transpose(Y_) * s; 
            
            // moving things one row up and deleting diagonal elements
            for(auto i =0; i < m_-1; i++)
            {
                for(auto j =0; j < i; j++)
                {
                    L_dots_[i][j] = L_dots_[i+1][j]; 
                }
                d_elements_[i] = d_elements_[i+1]; 
            }

            {
                Read<Vector> rw(dots); 
                auto r = range(dots); 
                
                for(auto i=r.begin(); i != r.end(); ++i)
                {
                    if(i==r.end()-1)
                        d_elements_[m_-1] = dots.get(i); 
                    else
                        L_dots_[m_-1][i] = dots.get(i); 
                }
            }
        }


        void buildM()
        {           
            // Doesn't work with dense matrices 
            DenseMatrix SS = theta_ * transpose(S_) * S_; 

            {
                Write<DenseMatrix> wM(M_); 
                Read<DenseMatrix> rS(SS); 

                auto r = row_range(M_); 
                auto c = col_range(M_); 

                // auto rowSS = row_range(SS); 
                // auto colSS = col_range(SS);  // here, we assume that SS has same distribution as M_ 

                for(auto i=r.begin(); i!=r.end(); ++i)
                {
                    for(auto j=c.begin(); j!=c.end(); ++j)
                    {
                        Scalar value = 0.0; 
                        
                        if(i < m_ && j < m_ && i==j) // top left block
                        {
                            value =  -1.0 * d_elements_[j]; 
                        }
                        else if(i >= m_ && j < m_ ) // bottom left  block 
                        {
                            value = L_dots_[i -m_][j]; 
                        }
                        else if(i < m_ && j >= m_ ) // top right  block 
                        {
                            value = L_dots_[j-m_][i]; 
                        }
                        else if(i >= m_ && j >= m_ ) // bottom right  block 
                        {
                            value = SS.get(i - m_, j - m_); 
                        }

                        M_.set(i,j, value); 
                    }
                }
            }


            // // because before, matrices are singular
            // if(current_m_ > m_)
            // {
            //     std::cout<<"------ in ------\n"; 
            //     DenseMatrix M_inv = local_values(local_size(M_).get(0), local_size(M_).get(1), 99.0); 

            //     if(IterativeSolver<Matrix, Vector> * ls = dynamic_cast<IterativeSolver<Matrix, Vector>*>(linear_solver_.get()))
            //         ls->atol(1e-15); 

            //     MatLinearSolver<DenseMatrix, DenseMatrix, Vector> mat_solver(linear_solver_); 
            //     mat_solver.get_inverse(M_, M_inv); 

            //     M_ = M_inv; 
            //     std::cout<<"------ out ------\n"; 
            // }

        }

    public:        

        bool constrained_solve(const Vector & x, const Vector & g, const Vector & lb, const Vector & ub, Vector & s) const override
        {
            Vector x_cp, c; 

            this->computeCauchyPoint(x, g, lb, ub, x_cp, c);
            this->compute_reduced_Newton_dir(x, x_cp, c, g, lb, ub, s); 

            return true; 
        }



    Scalar get_gb(const Vector & g, const SizeType & index) const
    {
        Vector g_local = local_values(1, 0); 

        {   // begin lock
            Write<Vector>  wd(g_local);
            Read<Vector> r1(g); 

            Range rr = range(g);

            Range rr_g_local = range(g_local);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(i==index)
                    g_local.set(rr_g_local.begin(), g.get(i)); 
            }
   
        } // end of lock

        return sum(g_local); 
    }



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        void computeCauchyPoint(const Vector &x, const Vector & g, const Vector & lb, 
                                const Vector & ub, Vector & x_cp, Vector & c) const
        {
            x_cp = x; 
            Scalar t_old, t, dt_min, dt, f_p, f_pp, z_b, g_b; 
            SizeType b; 
            Vector wbT; 

            const auto inf = std::numeric_limits<Scalar>::infinity(); 

            Vector break_points, feasible_set; // indexing and distribution as always
            cp_.compute_breakpoints(g, x, lb, ub, break_points); 

            Vector d; 
            cp_.get_d_corresponding_to_ti(break_points, g, d, 0.0); 
            cp_.get_initial_feasible_set(break_points, feasible_set); 

            // disp(feasible_set); 

            Vector sorted_break_points; 
            vec_unique_sort_serial(break_points, sorted_break_points, cp_.get_memory_size()); 

            Vector p; 
            // this two lines seem usless
            p = transpose(W_) * d; 
            c = local_values(local_size(W_).get(1), 1, 0); 

            f_p = -1.0 * dot(d, d); 
            f_pp = - theta_ * f_p - dot(p, M_ * p); 

            dt_min = -f_p/f_pp; 
            t_old = 0.0; 

            SizeType it_sorted = 1; 

            t = cp_.get_next_t(sorted_break_points, 0); 

            if(t==0)
            {
                t = cp_.get_next_t(sorted_break_points, 0); 
                it_sorted = 2;                    
            }

            dt = t - t_old; 


            bool repeated_index = cp_.get_global_active_index(break_points, feasible_set, t, b); 
            SizeType num_break_points = cp_.get_number_of_sorted_break_points(sorted_break_points); 
                

            // check logic here...
            while(dt_min >= dt && it_sorted < num_break_points)
            { 
                z_b = cp_.project_direction_on_boundary(x, d, ub, lb, x_cp, b); 
                g_b = get_gb(g, b); 

                c = c + dt * p; 

                DenseMatrix W_transpose_ = transpose(W_); 
                mat_get_col(W_transpose_, wbT, b); 

                Vector Mc, MwbT, Mp; 
                // if(current_m_ <= m_)
                // {
                    this->apply_M(c, Mc); 
                    this->apply_M(p, Mp);                 
                    this->apply_M(wbT, MwbT); 
                // }
                // else
                // {
                //     Mc = M_*c; 
                //     Mp = M_*p; 
                //     MwbT = M_*wbT; 
                // }

                f_p += (dt * f_pp) +  (g_b*g_b) + (theta_*g_b  * z_b); 
                f_p += g_b * dot(wbT, Mc); 

                f_pp -= (theta_ * (g_b*g_b)) + (2. * g_b * dot(wbT, Mp)) - ((g_b*g_b) * dot(wbT, MwbT)); 

                // TODO:: add checks 
                if(f_pp == 0 || !std::isfinite(f_pp))
                    break;   

                p += g_b * wbT; 
                cp_.zero_dir_component(d, b); 

                dt_min  = -f_p/f_pp;
                t_old   = t;

                // lets see 
                if(repeated_index)
                {
                    repeated_index = cp_.get_global_active_index(break_points, feasible_set, t, b); 
                    dt = t - t_old;
                }
                else
                {   
                    t = cp_.get_next_t(sorted_break_points, it_sorted);

                    if(t==inf || !std::isfinite(t))
                        break; 

                    repeated_index = cp_.get_global_active_index(break_points, feasible_set, t, b); 
                    it_sorted++; 

                    // TODO:: do checks for infinity.... 
                    dt = t - t_old;
                }

            }

            dt_min = std::max(dt_min, 0.0);
            t_old  = t_old + dt_min;

            cp_.add_d_to_x(x, x_cp, feasible_set, d, t_old); 

            c = c + dt_min*p;
        }



    // TODO:: return correction
    // returns true, if there is any free variable, otherwise return false
    // TODO:: investigate if you need new feasible set         
    bool compute_reduced_Newton_dir(const Vector & x,     const Vector & x_cp, const Vector & c, const Vector &g, 
                                    const Vector & lb,  const Vector & ub,  Vector & correction) const
    {
        Vector feasible_set; 

        correction = x_cp - x; 
        Vector help_g; 
        this->apply_H(correction, help_g); 
        Vector grad_quad_fun = -1.0 * (g + help_g); 

        // building feasible set 
        reduced_primal_method_.build_feasible_set(x_cp, ub, lb, feasible_set); 
        SizeType feasible_variables = sum(feasible_set); 

        Vector  s; 
    
        if(feasible_variables == 0) // all variables are feasible => perform Newton step on whole matrix
            return false; 
        else if(size(feasible_set).get(0)==feasible_variables)
        {
            this->apply_Hinv(grad_quad_fun, s); 
            Scalar alpha_star = reduced_primal_method_.compute_alpha_star(x_cp, lb, ub, s, feasible_set); 
            correction = correction + (alpha_star * s);    
            
        }
        else
        {
            this->apply_reduced_Hinv(feasible_set, grad_quad_fun, s); 

            Scalar alpha_star = reduced_primal_method_.compute_alpha_star(x_cp, lb, ub, s, feasible_set); 
        
            Vector s_prolongated; 
            reduced_primal_method_.prolongate_reduced_corr(s, feasible_set,  s_prolongated); 

            // final correction - both CP and Newton step 
            correction = correction + (alpha_star * s_prolongated); 
        }

        return true; 
    }







    protected: 
        // Sherman - Morison - Woodbury formula for computing inverse
        // this formula could be spped-up by using M_inv instead of M 
        virtual void apply_inverse_to_vec(const Vector &  g, const DenseMatrix & W, Vector & s) const 
        {
            Vector WTg = -1.0 * transpose(W) * g; 
            Vector v; 

            // if(current_m_ <= m_)
            // {
                this->apply_M(WTg, v); 
            // }
            // else
            // {
            //     v = M_*WTg; 
            // }

            DenseMatrix N  = 1.0/theta_ * transpose(W_) * W_; 
            
            DenseMatrix MN; 
            // if(current_m_ <= m_)
            // {
                this->apply_M(N, MN); 
            // }
            // else
            // {   
            //     MN = M_*N; 
            // }

            N = DenseMatrix(local_identity(local_size(MN))) - MN; 

            Vector N_inv_v = local_zeros(local_size(v)); 
            linear_solver_->solve(N, v, N_inv_v); 
            
            s = (1.0/theta_ )* g; 
            s -= 1.0/(theta_*theta_) * (W* N_inv_v); 
        }


        void apply_M(const Vector & v, Vector & result) const 
        {
            if(IterativeSolver<Matrix, Vector> * ls = dynamic_cast<IterativeSolver<Matrix, Vector>*>(linear_solver_.get()))
                    ls->atol(1e-15); 

            linear_solver_->solve(M_, v, result);
        }


        void apply_M(const DenseMatrix & RHS, DenseMatrix & result) const 
        {
            result = local_values(local_size(RHS).get(0), local_size(RHS).get(1), 0.0); 
            
            if(IterativeSolver<Matrix, Vector> * ls = dynamic_cast<IterativeSolver<Matrix, Vector>*>(linear_solver_.get()))
                    ls->atol(1e-15); 

            MatLinearSolver<DenseMatrix, DenseMatrix, Vector> mat_solver(linear_solver_); 
            mat_solver.solve(M_, RHS, result); 
        }



    private:
        // static_assert(!utopia::is_sparse<Matrix>::value, "LBFGS does not support sparse matrices.");

        SizeType m_; // memory size

        std::vector<Scalar> d_elements_; 
        std::vector<std::vector<Scalar> > L_dots_; 


        SizeType current_m_; // iteration number 
        Matrix H_;  // identity

        Scalar theta_; // scaling param for H0

        DenseMatrix Y_;  // matrix n \times m of grad diffs, Y_k = [y_{k-m}, ..., y_{k-1}], where y_k = g_{k+1} - g_k
        DenseMatrix S_;  // matrix n \times m of directions, S_k = [s_{k-m}, ..., s_{k-1}], where s_k = x_{k+1} - x_k
        
        DenseMatrix W_;  // W_k = [Y_k \theta S_k]
        DenseMatrix M_;  /* M_k  =   | -D       L_k^T  |
                                | L_k      \theta S_k^T S_k | */   // in contrast to the paper, this is not the inverse

        std::shared_ptr<LinSolver> linear_solver_;     /*!< Linear solver parameters. - Ideally LU  */  


        GeneralizedCauchyPoint<Matrix, Vector> cp_; 
        ReducedPrimalMethod<Matrix, Vector> reduced_primal_method_; 
};



}

#endif //UTOPIA_LBFGSB_HPP

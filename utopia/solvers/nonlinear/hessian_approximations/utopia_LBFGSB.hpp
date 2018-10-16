#ifndef UTOPIA_LBFGSB_HPP
#define UTOPIA_LBFGSB_HPP

#include "utopia_Core.hpp"

namespace utopia
{

template<class Matrix, class Vector>
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
            H0_ = local_identity(n, n);

            if(cp_.get_memory_size()==-1)
                cp_.set_memory_size(size(x).get(0));

            this->initialized(true);
            current_m_ = 0; 

            // TODO:: recheck 
            W_ = values(size(x).get(0), 1, 0.0); 
            M_ = zeros(1, 1); 

            // TODO:: investigate proper value of theta
            theta_ = 1.0; 

            return true;
        }


        /**
         * @brief      Changes linear solver used inside of nonlinear-solver. 
         *
         * @param[in]  linear_solver  The linear solver
         */
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

            // TODO:: check if update is SPD ... 

            Scalar denom    = dot(y,s); 
            Scalar nom      = dot(y,y);

            // if denom > eps, hessian approx. should be positive semidefinite
            if(denom < 1e-12)
            {
                utopia_warning("L-BFGS-B: Curvature condition not satified. Skipping update. ")
                return false; 
            }

            theta_ = nom/denom; 

            // this needs to be done durign first it, in order to initialize matrices 
            if(current_m_ == 0)
            {
                this->init_mat_from_vec(Y_, y); 
                this->init_mat_from_vec(S_, s); 
            }
            // here we are just building first m elements of matrices
            else if(current_m_ < m_)
            {
                this->add_col_to_mat(Y_, y); 
                this->add_col_to_mat(S_, s); 
            }
            // here, we are removing first col, shifting, appending new 
            else
            {  
                this->shift_cols_left_replace_last_col(Y_, y); 
                this->shift_cols_left_replace_last_col(S_, s); 

            }

            this->buildW(); 
            this->buildM(); 

            current_m_++; 

            return true;
        }

        virtual bool apply_Hinv(const Vector & g, Vector & s) const override
        {
            std::cerr << "--- not implemented yet---- \n";
            return true;
        }


        virtual bool apply_H(const Vector & v, Vector & result) const override
        {
            Vector Y_v = transpose(Y_) * v;
            Vector S_v = theta_ * transpose(S_) * v;

            Vector p =  Vector(Blocks<Vector>(
            {
                make_ref(Y_v), make_ref(S_v)
            }));

            this->apply_M(p, result); 
            result = (theta_ * v) - W_*result; 

            return false; 
        }


        virtual Matrix & get_Hessian() override
        {
            
            utopia_warning("LBFGS::get_Hessian returns dense matrix ...."); 
            return H0_;
        }

        virtual Matrix & get_Hessian_inv() override
        {
            std::cerr << "--- not implemented yet---- \n";
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


    private:
        void add_col_to_mat(Matrix & M, const Vector & col) const
        {
            const Matrix M_old = M;
            M = values(size(M_old).get(0), size(M_old).get(1) + 1, 0.0);

            {   // begin lock

                Write<Matrix>  write(M);

                Read<Matrix>   read(M_old);
                Read<Vector>      read_vec(col);

                Range r = row_range(M_old);
                Range c = col_range(M_old);

                Range r_new = row_range(M);
                Range c_new = col_range(M);

                for (SizeType i = r.begin(); i != r.end(); ++i) {
                    for (SizeType j = c.begin(); j != c.end(); ++j)
                        M.set(i, j, M_old.get(i, j));
                }

                for (SizeType i = r_new.begin(); i != r_new.end(); ++i)
                    M.set(i, c_new.end() - 1, col.get(i));

            } // end of lock
        }

        void shift_cols_left_replace_last_col(Matrix & M, const Vector & col) const
        {
            const Matrix M_old = M;

            {   // begin lock

                Write<Matrix>  write(M);

                Read<Matrix>   read(M_old);
                Read<Vector>      read_vec(col);

                Range r = row_range(M_old);
                Range c = col_range(M_old);

                for (SizeType i = r.begin(); i != r.end(); ++i) {
                    for (SizeType j = c.begin(); j != c.end() - 1; ++j) {
                        M.set(i, j, M_old.get(i, j + 1));
                    }
                    M.set(i, c.end() - 1, col.get(i));
                }

            } // end of lock
        }

        void init_mat_from_vec(Matrix & M, const Vector & col) const 
        {
            M = values(size(col).get(0), 1, 0.0);

            {   // begin lock

                Write<Matrix>  write(M);
                Read<Vector>   read_vec(col);

                Range r = row_range(M);
                Range c = col_range(M);

                for (SizeType i = r.begin(); i != r.end(); ++i)
                    M.set(i, c.end()-1, col.get(i));

            } // end of lock
        }

        void buildW()
        {
            Matrix S_theta = theta_ * S_; 

            W_ = Matrix(Blocks<Matrix>(1, 2,
            {
                make_ref(Y_), make_ref(S_theta)
            }));
        }


        void buildM()
        {
            SizeType n = (current_m_ < m_)? current_m_ + 1: m_; 

            Vector d = zeros(n); 

            if(size(S_)!=size(Y_) || n != size(Y_).get(1)){
                std::cout<<"LBFGSB::buildM:: sizes do not match .... \n"; 
                return; 
            }

            std::vector<Scalar> diag_dots(n); 

            for(auto i=0; i < n; i++)
            {
                // Vector v1 = local_zeros(local_size(S_)), v2 = local_zeros(local_size(S_)); 
                Vector v1, v2; 
                mat_get_col(S_, v1, i);
                mat_get_col(Y_, v2, i);
                diag_dots[i] = -1.0 * dot(v1, v2); 
            }

            {   // begin lock

                Write<Vector>   read_vec(d);
                Range r = range(d);

                for (SizeType i = r.begin(); i != r.end(); ++i)
                        d.set(i, -1.0 * diag_dots[i]);

            } // end of lock

            Matrix D = diag(d); 


            std::vector<std::vector<Scalar> > matrixL(n, std::vector<Scalar>(n));

            for(auto i=0; i < n; i++)
            {
                for(auto j=0; j < n; j++)
                {
                    if(i > j)
                    {
                        Vector v1, v2; 
                        mat_get_col(S_, v1, i);
                        mat_get_col(Y_, v2, j);
                        matrixL[i][j] = dot(v1, v2); 
                    }
                }
            }

            Matrix L = zeros(n, n); 

            {   // begin lock
                Write<Matrix>   write_mat(L);

                Range rr = row_range(L);
                Range cc = col_range(L);

                for (SizeType r = rr.begin(); r != rr.end(); ++r)
                {
                    for (SizeType c = cc.begin(); c != cc.end(); ++c)
                    {
                        if(r > c)
                            L.set(r, c,  matrixL[r][c]);
                    }
                }

            } // end of lock

            Matrix SS = theta_ * transpose(S_) * S_; 
            Matrix LT = transpose(L); 


            M_ = Matrix(Blocks<Matrix>(2, 2,
            {
                make_ref(D), make_ref(LT), 
                make_ref(L), make_ref(SS), 
            }));

        }

        void apply_M(const Vector & v, Vector & result) const 
        {
            linear_solver_->solve(M_, v, result);
        }


        void apply_M(const Matrix & RHS, Matrix & result) const 
        {
            // TO BE DONE:: 
            // linear_solver_->solve(M_, v, result);
            // MatLinearSolver<>
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

                Matrix W_transpose_ = transpose(W_); 
                mat_get_col(W_transpose_, wbT, b); 

                Vector Mc, MwbT, Mp; 
                this->apply_M(c, Mc); 
                this->apply_M(p, Mp);                 
                this->apply_M(wbT, MwbT); 

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
        Matrix W_reduced; 
        Vector g_reduced;
        Vector feasible_set; 

        Vector Mc; 
        this->apply_M(c, Mc); 

        Vector global_grad   = g + (theta_*(x_cp-x)) - (W_*(Mc)); 
        Matrix W_T = transpose(W_);  

        reduced_primal_method_.build_reduced_quantities(lb, ub, x_cp, global_grad, W_T, feasible_set,  W_reduced, g_reduced); 


        if(sum(feasible_set)==0)
        {
            correction = x_cp - x; 
            return false; 
        }

        Vector WR_gR = W_reduced * g_reduced; 

        Vector v; 
        this->apply_M(WR_gR, v); 

        Matrix N = W_T * transpose(W_T); 
        N  = 1/theta_ * N; 
        N = M_ * N; 
        N = Matrix(local_identity(local_size(N))) - N; 

        Vector s = local_zeros(local_size(v)); 
        linear_solver_->solve(N, v, s); 
        

        Vector  du = (-1.0/theta_ )* g_reduced; 
                du -= 1.0/(theta_*theta_) * transpose(W_reduced) * s; 

        Scalar alpha_star = reduced_primal_method_.compute_alpha_star(x_cp, lb, ub, du, feasible_set); 
        du *= alpha_star; 

        Vector x_bar = x_cp; 
        Vector du_prolongated; 

        reduced_primal_method_.prolongate_reduced_corr(du, feasible_set,  du_prolongated); 

        x_bar = x_bar + du_prolongated; 
        correction = x_bar - x; 

        return true; 
    }



    private:
        static_assert(utopia::is_sparse<Matrix>::value, "LBFGS does not support dense matrices.");

        SizeType m_; // memory size

        SizeType current_m_; // iteration number 
        Matrix H0_;  // identity

        Scalar theta_; // scaling param for H0

        Matrix Y_;  // matrix n \times m of grad diffs, Y_k = [y_{k-m}, ..., y_{k-1}], where y_k = g_{k+1} - g_k
        Matrix S_;  // matrix n \times m of directions, S_k = [s_{k-m}, ..., s_{k-1}], where s_k = x_{k+1} - x_k
        
        Matrix W_;  // W_k = [Y_k \theta S_k]
        Matrix M_;  /* M_k  =   | -D       L_k^T  |
                                | L_k      \theta S_k^T S_k | */   // in contrast to the paper, this is not the inverse

        std::shared_ptr<LinSolver> linear_solver_;     /*!< Linear solver parameters. - Ideally LU  */  


        GeneralizedCauchyPoint<Matrix, Vector> cp_; 
        ReducedPrimalMethod<Matrix, Vector> reduced_primal_method_; 
};



}

#endif //UTOPIA_LBFGSB_HPP

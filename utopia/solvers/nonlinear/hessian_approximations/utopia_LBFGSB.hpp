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
    typedef utopia::LSStrategy<Matrix, Vector>      LSStrategy; 

    public:

        LBFGSB( const SizeType & m, 
                const std::shared_ptr <LinSolver> &linear_solver = std::make_shared<ConjugateGradient<Matrix, Vector> >(),
                const std::shared_ptr <LSStrategy> &line_search = std::make_shared<SimpleBacktracking<Matrix, Vector> >()):
                m_(m), cp_memory_(-1), current_m_(0), linear_solver_(linear_solver), line_search_strategy_(line_search)
        {

        }

        virtual bool initialize(Function<Matrix, Vector> &fun, const Vector &x) override
        {
            SizeType n = local_size(x).get(0);
            H0_ = local_identity(n, n);

            if(cp_memory_==-1)
                cp_memory_ = size(x).get(0);

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
         * @brief      Sets strategy for computing step-size. 
         *
         * @param[in]  strategy  The line-search strategy.
         *
         * @return     
         */
        virtual bool set_line_search_strategy(const std::shared_ptr<LSStrategy> &strategy)
        {
          line_search_strategy_ = strategy; 
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

            theta_ = dot(y,y)/dot(y,s); 

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

        virtual bool apply_Hinv(const Vector & /* g */, Vector & /*s */) const override
        {
            // TODO
            return true;
        }


        virtual Matrix & get_Hessian() override
        {
            std::cerr << "--- not implemented yet---- \n";
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

        void set_Cauchy_point_memory_size(const SizeType & m)
        {
            cp_memory_ = m;
        }

        SizeType get_Cauchy_point_memory_size() const 
        {
            return cp_memory_;
        }



    private:
        /**
         * @brief Appends vector into matrix as new last column
         * 
         * @param M Matrix to be appended to
         * @param col Vector to be added
         */
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

        /**
         * @brief Add shift columns of the matrix to the left by one col. 
         * Then, replace last column with new vector
         * 
         * @param M Matrix to be shifted
         * @param col Vector to be added to the last col
         */
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

        /**
         * @brief Creates matrix from the vector. Basically new matrix has dimension n\times 1 
         * 
         * @param M matrix
         * @param col vecotr
         */
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

        /**
         * @brief Build matrix W
         */
        void buildW()
        {
            Matrix S_theta = theta_ * S_; 

            W_ = Matrix(Blocks<Matrix>(1, 2,
            {
                make_ref(Y_), make_ref(S_theta)
            }));
        }


        /**
         * @brief Build matrix M from vectors Y and S
         * 
         */
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

        void apply_M(const Vector & v, Vector & result)
        {
            linear_solver_->solve(M_, v, result);
        }


    public:        

        /**
         * @brief Computes breakpoints 
         * @details TODO:: add TR radius 
         * 
         * @param g gradient
         * @param x current iterate
         * @param lb lower bound
         * @param ub uppre bound
         * @param t breakpoints - distributed vector for each point 
         */
        void compute_breakpoints(const Vector & g, const Vector & x, const Vector & lb, const Vector & ub, Vector &t)
        {
            auto inf = std::numeric_limits<Scalar>::infinity(); 

            if(empty(t) || local_size(t)!=local_size(x))
                t = local_values(local_size(x).get(0), inf);

            // TODO:: add checks if there are not all bounds available 

            {
              Read<Vector> r_ub(ub), r_lb(lb), r_x(x), r_d(g);
              Write<Vector> wt(t); 

              each_write(t, [ub, lb, x, g, inf](const SizeType i) -> double { 
                          Scalar li =  lb.get(i); Scalar ui =  ub.get(i); Scalar xi =  x.get(i);  Scalar gi =  g.get(i);  
                          if(gi < 0)
                            return (xi - ui)/gi; 
                          else if(gi > 0)
                            return (xi - li)/gi; 
                        else 
                            return inf; 
            }  );
          }
        }



        /**
         * @brief Determines grad. descent direction, such that lb and ub are not violated
         * 
         * @param t break points
         * @param g gradient 
         * @param d search direction
         * @param t_current breakpoint to be found
         */
        void get_d_corresponding_to_ti(const Vector & t, const Vector & g, Vector &d, const Scalar & t_current)
        {
            d = -1.0 * g; 

            {   // begin lock
                Write<Vector>  wd(d);
                Read<Vector>   rt(t);

                Range rr = range(d);

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    Scalar ti = t.get(i); 
                    if(std::abs(ti - t_current) < 1e-12)
                        d.set(i, 0.0); 
                }

            } // end of lock
        }

    void get_initial_feasible_set(const Vector & break_points, Vector & feasible_set)
    {
        feasible_set = local_values(local_size(break_points).get(0), 0.0); 

        {   // begin lock
            Read<Vector>  wd(break_points);
            Write<Vector> r(feasible_set); 

            Range rr = range(break_points);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(break_points.get(i) > 0)
                    feasible_set.set(i, 1); 
            }
   
        } // end of lock
    }



    bool get_global_active_index(const Vector & break_points, Vector & feasible_set, const Scalar & t_current, SizeType & index)
    {
        SizeType prev_index = index; 
        SizeType counter_global_sum;  
        Scalar counter=0.0; 

        // TODO:: put inf
        Vector indices = local_values(1, 9e9); 
        Vector counter_vec = local_values(1, 0); 

        {   // begin lock
            Read<Vector>  wd(break_points);
            Read<Vector>  rv(feasible_set); 

            Range rr = range(break_points);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                Scalar value = break_points.get(i); 
                if(std::abs(value - t_current) < 1e-12 && feasible_set.get(i)==1)
                {
                    indices.set(0,i); 
                    break; 
                }
            }

            index = min(indices); 

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(i==index)
                    feasible_set.set(i, 0); 
            }            


            // horrible solution - looops should be merged 
            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                Scalar value = break_points.get(i); 
                if(std::abs(value - t_current) < 1e-12 && feasible_set.get(i)==1)
                {
                    counter++; 
                }
            }    


            Write<Vector> wvv(counter_vec); 

            Range r_counter = range(counter_vec);
            for (SizeType i = r_counter.begin(); i != r_counter.end(); ++i)
                counter_vec.set(i, counter); 


        } // end of lock

        counter_global_sum = sum(counter_vec); 

        std::cout<<"counter_global_sum  \n"; 

        return (counter_global_sum>1) ? true: false; 
    }




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        void computeCauchyPoint(const Vector &x, const Vector & g, const Vector & lb, const Vector & ub, Vector & x_cp, Vector & c)
        {
            x_cp = x; 
            Scalar t_old, t, dt_min, dt, f_p, f_pp; 
            SizeType b; 

            Vector break_points, feasible_set; // indexing and distribution as always
            this->compute_breakpoints(g, x, lb, ub, break_points); 

            Vector d; 
            this->get_d_corresponding_to_ti(break_points, g, d, 0.0); 
            this->get_initial_feasible_set(break_points, feasible_set); 

            // disp(feasible_set); 

            Vector sorted_break_points; 
            vec_unique_sort_serial(break_points, sorted_break_points, cp_memory_); 

            Vector p; 
            // this two lines seem usless
            p = transpose(W_) * d; 
            c = zeros(local_size(W_).get(1), 1); 

            f_p = -1.0 * dot(d, d); 
            f_pp = - theta_ * f_p - dot(p, M_ * p); 

            dt_min = -f_p/f_pp; 
            t_old = 0.0; 

            Vector t_help = local_values(1, 0.0); 

            // this is horrible solution, but lets fix it later 
            {
                Read<Vector> rv(sorted_break_points); 
                auto rr = range(sorted_break_points);
                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {
                    if(i==0)
                        t_help.add(0, sorted_break_points.get(i)); 
                }
            }

            t = sum(t_help); 
            dt = t - t_old; 

            bool repeated_index = get_global_active_index(break_points, feasible_set, t, b); 











        }







    private:
        static_assert(utopia::is_sparse<Matrix>::value, "BFGS does not support sparse matrices.");

        SizeType m_; // memory size
        SizeType cp_memory_; // memory size

        SizeType current_m_; // iteration number 
        Matrix H0_;  // identity

        Scalar theta_; // scaling param for H0

        Matrix Y_;  // matrix n \times m of grad diffs, Y_k = [y_{k-m}, ..., y_{k-1}], where y_k = g_{k+1} - g_k
        Matrix S_;  // matrix n \times m of directions, S_k = [s_{k-m}, ..., s_{k-1}], where s_k = x_{k+1} - x_k
        
        Matrix W_;  // W_k = [Y_k \theta S_k]
        Matrix M_;  /* M_k  =   | -D       L_k^T  |
                                | L_k      \theta S_k^T S_k | */   // in contrast to the paper, this is not the inverse

        std::shared_ptr<LinSolver> linear_solver_;     /*!< Linear solver parameters. - Ideally LU  */  
        std::shared_ptr<LSStrategy> line_search_strategy_;     /*!< Strategy used in order to obtain step \f$ \alpha_k \f$ */  
};



}

#endif //UTOPIA_LBFGSB_HPP

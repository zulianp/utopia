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

            Scalar denom    = dot(y,s); 
            Scalar nom      = dot(y,y);

            // if(denom > 1e-12 * nom)
            // {
            //     utopia_warning("L-BFGS-B: Curvature condition not satified. Skipping update. ")
            //     return false; 
            // }

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


    public:        

        bool constrained_solve(const Vector & x, const Vector & g, const Vector & lb, const Vector & ub, Vector & s) const override
        {
            Vector x_cp, c; 

            this->computeCauchyPoint(x, g, lb, ub, x_cp, c);
            this->compute_reduced_Newton_dir(x, x_cp, c, g, lb, ub, s); 

            return true; 
        }




        void compute_breakpoints(const Vector & g, const Vector & x, const Vector & lb, const Vector & ub, Vector &t) const
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



        void get_d_corresponding_to_ti(const Vector & t, const Vector & g, Vector &d, const Scalar & t_current) const
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

    void get_initial_feasible_set(const Vector & break_points, Vector & feasible_set) const
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



    bool get_global_active_index(const Vector & break_points, Vector & feasible_set, const Scalar & t_current, SizeType & index) const
    {
        Scalar counter=0.0; 

        // TODO:: put inf
        Vector indices = local_values(1, 9e9); 
        Vector counter_vec = local_values(1, 0); 

        {   // begin lock
            Read<Vector>  wd(break_points);
            Read<Vector>  rv(feasible_set); 
            Write<Vector>  rv2(indices); 

            Range rr = range(break_points);
            Range ind_range = range(indices);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                Scalar value = break_points.get(i); 
                if(std::abs(value - t_current) < 1e-12 && feasible_set.get(i)==1)
                {
                    indices.set(ind_range.begin(), i); 
                    break; 
                }
            }
        }

        index = min(indices); 

        if(index<0 || index > size(break_points).get(0))
            utopia_error("L-BFGS-B::get_global_active_index: index not valid. "); 


        {
            Write<Vector>  rv(feasible_set); 
            Range rr = range(feasible_set);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(i==index)
                    feasible_set.set(i, 0); 
            }            

        }

        {
            Read<Vector>  rv(feasible_set); 
            Range rr = range(feasible_set);

            // horrible solution - looops should be merged 
            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                Scalar value = break_points.get(i); 
                if(std::abs(value - t_current) < 1e-12 && feasible_set.get(i)==1)
                {
                    counter++; 
                }
            }
        }    

        {
            Write<Vector> wvv(counter_vec); 
            Range r_counter = range(counter_vec);

            for (SizeType i = r_counter.begin(); i != r_counter.end(); ++i)
                counter_vec.set(i, counter); 

        } // end of lock

        SizeType counter_global_sum = sum(counter_vec); 
        return (counter_global_sum>1) ? true: false; 
    }




    Scalar project_direction_on_boundary(const Vector & x, const Vector & d, const Vector & ub, const Vector & lb, Vector & x_cp, const SizeType & active_index) const
    {
        Vector x_local = local_values(1, 0.0); 
        Scalar val=0; 

        {   // begin lock
            Write<Vector>  wd(x_cp);

            Read<Vector> r1(ub); 
            Read<Vector> r2(lb); 
            Read<Vector> r3(d); 
            Read<Vector> r4(x); 

            Range rr = range(x_cp);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(i==active_index)
                {
                    val = (d.get(i)>0)? ub.get(i) : lb.get(i); 
                    x_cp.set(i, val); 
                }
            }
   
        } // end of lock


        ///////////////////////////////////////////////////////////
        {   // begin lock
            Write<Vector>  wx(x_local);
            Range rr_x_local = range(x_local);

            x_local.set(rr_x_local.begin(), val); 
   
        } // end of lock

        return sum(x_local); 
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


    void zero_dir_component(Vector & d, const SizeType & index) const
    {
        {   // begin lock
            Write<Vector>  wd(d);
            Range rr = range(d);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(i==index)
                    d.set(i, 0); 
            }
   
        } // end of lock

    }


    // use approxeq for all stuff where u compare 
    void add_d_to_x(const Vector & x, Vector & x_cp, const Vector & feasible_set,  const Vector & d, const Scalar & tau) const
    {

        {   // begin lock
            Write<Vector>  wd(x_cp);
            Read<Vector>  r1(x);
            Read<Vector>  r2(feasible_set);
            Read<Vector>  r3(d);

            Range rr = range(x_cp);

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(feasible_set.get(i) == 1.0)
                {
                    x_cp.set(i, x.get(i) + tau * d.get(i)); 
                }
                    
            }
   
        } // end of lock

    }



    Scalar get_next_t(const Vector & sorted_break_points, const SizeType & index) const
    {
        Vector t_help = local_values(1, 0.0); 
        Scalar value=0.0; 

        // this is horrible solution, but lets fix it later 
        {
            Read<Vector> r(sorted_break_points); 

            auto rr = range(sorted_break_points);
            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                if(i==index)
                    value = sorted_break_points.get(i); 
            }
        }

        // this is horrible solution, but lets fix it later 
        {
            Write<Vector> w(t_help); 

            auto rr = range(t_help);
            for (SizeType i = rr.begin(); i != rr.end(); ++i)
                    t_help.set(i, value); 
        }        

        return sum(t_help); 
    }



    SizeType get_number_of_sorted_break_points(const Vector & sorted_break_points) const
    {
        Vector help = local_values(1, 0.0); 
        SizeType val = size(sorted_break_points).get(0); 

        // this is horrible solution, but lets fix it later 
        {
            Write<Vector> w(help); 

            if(mpi_world_rank()==0)
                help.set(0, val); 
        }

        return sum(help); 
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
            c = local_values(local_size(W_).get(1), 1, 0); 

            f_p = -1.0 * dot(d, d); 
            f_pp = - theta_ * f_p - dot(p, M_ * p); 

            dt_min = -f_p/f_pp; 
            t_old = 0.0; 

            SizeType it_sorted = 1; 

            t = get_next_t(sorted_break_points, 0); 

            if(t==0)
            {
                t = get_next_t(sorted_break_points, 0); 
                it_sorted = 2;                    
            }

            dt = t - t_old; 


            bool repeated_index = get_global_active_index(break_points, feasible_set, t, b); 
            SizeType num_break_points = get_number_of_sorted_break_points(sorted_break_points); 
                

            // check logic here...
            while(dt_min >= dt && it_sorted < num_break_points)
            { 
                z_b = project_direction_on_boundary(x, d, ub, lb, x_cp, b); 
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
                this->zero_dir_component(d, b); 

                dt_min  = -f_p/f_pp;
                t_old   = t;

                // lets see 
                if(repeated_index)
                {
                    repeated_index = get_global_active_index(break_points, feasible_set, t, b); 
                    dt = t - t_old;
                }
                else
                {   
                    t = get_next_t(sorted_break_points, it_sorted);

                    if(t==inf || !std::isfinite(t))
                        break; 

                    repeated_index = get_global_active_index(break_points, feasible_set, t, b); 
                    it_sorted++; 

                    // TODO:: do checks for infinity.... 
                    dt = t - t_old;
                }

            }

            dt_min = std::max(dt_min, 0.0);
            t_old  = t_old + dt_min;

            this->add_d_to_x(x, x_cp, feasible_set, d, t_old); 

            c = c + dt_min*p;
        }



        void build_reduced_quantities(  const Vector &lb, const Vector & ub, const Vector & x_cp,
                                        const Vector & g, const Matrix & H, 
                                        Vector & feasible_set,  Matrix & H_reduced, Vector & g_reduced) const
        {
            feasible_set = local_values(local_size(lb).get(0), 0.); 

            // std::cout<<"build_reduced_quantities: 1 ---- \n"; 

            SizeType local_feasible_set = 0; 

            {
                Write<Vector>  rv(feasible_set); 

                Read<Vector>  rv1(ub); 
                Read<Vector>  rv2(lb); 
                Read<Vector>  rv3(x_cp); 

                auto rr = range(feasible_set); 

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {   
                    // TODO:: put approx eq
                    if(lb.get(i) != x_cp.get(i) &&  ub.get(i) != x_cp.get(i)){
                        feasible_set.set(i, 1); 
                        local_feasible_set++; 
                    }
                }       
            }

            // std::cout<<"build_reduced_quantities: 2 ---- \n"; 

            // std::cout<<"----- feasible set:   \n"; 
            // disp(feasible_set); 

            // we can just get things directly 
            if(size(feasible_set).get(0)==sum(feasible_set))
            {
                g_reduced = g; 
                H_reduced = H; 

                return; 
            }


            g_reduced = local_zeros(local_feasible_set); 

            {
                SizeType local_counter = 0; 

                Read<Vector>  rf(feasible_set); 
                Read<Vector>  rg(g); 

                Write<Vector>  wg(g_reduced); 
                auto rr = range(feasible_set); 
                auto range_reduced = range(g_reduced); 

                for (SizeType i = rr.begin(); i != rr.end(); ++i)
                {   
                    // TODO:: put approx eq
                    if(feasible_set.get(i)==1)
                    {
                        g_reduced.set(range_reduced.begin() + local_counter, g.get(i)); 
                        local_counter++;
                    }
                }       
            }

            // std::cout<<"build_reduced_quantities: 3 ---- \n"; 

            // TODO:: use dense matrices for this kind of things 
            H_reduced  = local_values(local_size(H).get(0), local_feasible_set, 0.0); 

            // std::cout<<"-------- H -------- \n"; 
            // disp(H); 

            // std::cout<<"H_reduced.local(0): "<< local_size(H_reduced).get(0) << "  local_size(H_reduced).get(1):  "<< local_size(H_reduced).get(1) << "  \n"; 
            // std::cout<<"feasible_set: "<< local_size(feasible_set).get(0) << "  \n"; 
            // std::cout<<"H.local(0): "<< local_size(H).get(0) << "  local_size(H).get(1):  "<< local_size(H).get(1) << "  \n"; 


            if(local_size(feasible_set).get(0) != local_size(H).get(1))
                utopia_error("feasible_set, H: local sizes do not match .... \n"); 

            if(local_size(H_reduced).get(0) != local_size(H).get(0))
                utopia_error("H_reduced, H: local sizes do not match .... \n");             

            {
                Read<Vector>  rf(feasible_set); 
                Read<Matrix>  rM(H); 

                Write<Matrix>  wg(H_reduced); 

                auto row_original = row_range(H); 
                auto col_original = col_range(H); 

                // auto row_reduced = row_range(H_reduced); 
                auto col_reduced = col_range(H_reduced); 

                SizeType local_counter = 0; 

                // for (SizeType c = col_original.begin(); c != col_original.end(); ++c)
                for (SizeType r = row_original.begin(); r != row_original.end(); ++r)
                {   

                    for (SizeType c = col_original.begin(); c != col_original.end(); ++c)
                    {            
                        if(c==col_original.begin())
                            local_counter=0;

                        // std::cout<<"r: "<< r << "  c: "<< c << "  \n"; 

                        // TODO:: put approx eq
                        if(feasible_set.get(c)==1)
                        {
                            H_reduced.set(r, col_reduced.begin() + local_counter,  H.get(r, c)); 
                            local_counter++;
                        }
                    }    


                }  

            }

            // std::cout<<"-------- H_reduced -------- \n"; 
            // disp(H_reduced); 

        }


    Scalar compute_alpha_star(const Vector & x_cp, const Vector & lb, const Vector & ub, const Vector & d, const Vector & feasible_set) const
    {
        Vector alpha_stars = local_values(local_size(feasible_set).get(0), 1.0);  

        {
            Read<Vector>  rf(feasible_set); 
            Read<Vector>  rd(d); 
            Read<Vector>  ru(ub); 
            Read<Vector>  rlb(lb); 

            Write<Vector> wv(alpha_stars); 

            auto rr = range(alpha_stars); 

            for (SizeType i = rr.begin(); i != rr.end(); ++i)
            {
                // TODO:: put approx eq
                if(feasible_set.get(i)==1)
                {
                    Scalar val = 1.0; 
                    Scalar di = d.get(i);  Scalar xi = x_cp.get(i); Scalar li = lb.get(i);  Scalar ui = ub.get(i); 

                    if(di!=0)
                        val = (di>0) ? (ui - xi)/di : (li - xi)/di; 

                    // checks for nans and infs
                    if(val < 1.0 && std::isfinite(val) && val > -1e9)
                        alpha_stars.set(i, val); 

                }
            }     
        }

        return min(alpha_stars); 
    }




    void prolongate_reduced_corr(const Vector & x_reduced,  const Vector & feasible_set,  Vector & x_prolongated) const
    {
        x_prolongated = local_values(local_size(feasible_set).get(0), 0.0); 


        {
            Write<Vector>  wx(x_prolongated); 

            Read<Vector>   rv(x_reduced); 
            Read<Vector>   rf(feasible_set); 

            auto rr = range(feasible_set); 
            auto rr_red = range(x_reduced); 

            SizeType counter = 0; 

            for(auto i = rr.begin(); i != rr.end(); ++i)
            {
                if(feasible_set.get(i)==1)
                {
                    Scalar val = x_reduced.get(rr_red.begin() + counter); 
                    x_prolongated.set(i, val); 
                    counter++; 
                }
            }
        }
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

        this->build_reduced_quantities(lb, ub, x_cp, global_grad, W_T, feasible_set,  W_reduced, g_reduced); 


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


        Scalar alpha_star = this->compute_alpha_star(x_cp, lb, ub, du, feasible_set); 
        du *= alpha_star; 

        Vector x_bar = x_cp; 
        Vector du_prolongated; 

        this->prolongate_reduced_corr(du, feasible_set,  du_prolongated); 

        x_bar = x_bar + du_prolongated; 
        correction = x_bar - x; 

        return true; 
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

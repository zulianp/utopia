/*
* @Author: alenakopanicakova
* @Date:   2016-03-29
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/

#ifndef UTOPIA_MULTIGRID_HPP
#define UTOPIA_MULTIGRID_HPP
#include "utopia_Smoother.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_Utils.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_Level.hpp"
#include "utopia_MultiLevelBase.hpp"
#include <ctime>


namespace utopia 
{
    /**
     * @brief      The class for Linear Multigrid solver. 
     *
     * @tparam     Matrix  
     * @tparam     Vector  
     */
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Multigrid : public MultiLevelBase<Matrix, Vector>, 
                      public IterativeSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::LinearSolver<Matrix, Vector>        Solver;
        typedef utopia::IterativeSolver<Matrix, Vector>     IterativeSolver;
        typedef utopia::Smoother<Matrix, Vector>            Smoother;
        typedef utopia::Level<Matrix, Vector>               Level;
        typedef utopia::Transfer<Matrix, Vector>            Transfer;

        typedef struct 
        {
            std::vector<Vector> r_h, r_H, c_H, c_h; 

            void init(const int num_levels)
            {
                r_h.resize(num_levels);
                r_H.resize(num_levels);
                c_H.resize(num_levels);
                c_h.resize(num_levels);
            }

            inline bool valid(const std::size_t num_levels) const {
                return r_h.size() == num_levels;
            }

        } LevelMemory;

        LevelMemory memory;


    // #define BENCHMARKING_mode
   // #define CHECK_NUM_PRECISION_mode

    public:
       static const int V_CYCLE = 1;
       static const int W_CYCLE = 2;

       /**
        * @brief      Multigrid class
        *
        * @param[in]  smoother       The smoother.
        * @param[in]  direct_solver  The direct solver for coarse level. 
        */
        Multigrid(const std::shared_ptr<Smoother> &smoother = std::shared_ptr<Smoother>(), 
                  const std::shared_ptr<Solver> &direct_solver = std::shared_ptr<Solver>(),
                  const Parameters params = Parameters())
        : _smoother(smoother), _direct_solver(direct_solver), perform_galerkin_assembly_(true), use_line_search_(false)
        {
            set_parameters(params); 
        }

        virtual ~Multigrid(){} 
        
        void set_parameters(const Parameters params) override
        {
            IterativeSolver::set_parameters(params); 
            MultiLevelBase<Matrix, Vector>::set_parameters(params); 
            _smoother->set_parameters(params); 

        }

        /*! @brief if overriden the subclass has to also call this one first
         */
        virtual void update(const std::shared_ptr<const Matrix> &op) override
        {
            IterativeSolver::update(op);
            if(perform_galerkin_assembly_) {
                this->galerkin_assembly(op);
            }
        }


        virtual bool apply(const Vector &rhs, Vector & x_0) override
        {
            return solve(rhs, x_0);
        }

        /**
         * @brief      The solve function for multigrid method. 
         *
         * @param[in]  rhs   The right hand side.
         * @param      x_0   The initial guess. 
         *
         */
        virtual bool solve(const Vector &rhs, Vector &x_0)
        {
            Vector r, D_inv; 
            SizeType l = this->num_levels(); 

            memory.init(l);

            Scalar r_norm, r0_norm, rel_norm;
            SizeType it = 0; 

            bool converged = false; 

            r = rhs - level(l-1).A() * x_0; 

            r_norm = norm2(r); 
            r0_norm = r_norm; 


            this->init_solver("Multigrid", {" it. ", "|| r_N ||", "r_norm" }); 

            if(this->verbose())
                PrintInfo::print_iter_status(it, {r_norm, 1}); 
            it++; 

#ifdef BENCHMARKING_mode
                this->atol(1e-15); 
#endif

            if(this->cycle_type() == FULL_CYCLE)
                this->max_it(1); 

            while(!converged)
            {            
                if(this->cycle_type() == MULTIPLICATIVE_CYCLE)
                    multiplicative_cycle(rhs, l, x_0); 
                else if(this->cycle_type() == FULL_CYCLE)
                    full_cycle(rhs, l, x_0); 
                else
                    std::cout<<"ERROR::UTOPIA_MG<< unknown MG type... \n"; 


#ifdef CHECK_NUM_PRECISION_mode
                    if(has_nan_or_inf(x_0))
                    {
                        x_0 = local_zeros(local_size(x_0));
                        return true; 
                    }
#endif    

                r = rhs - level(l-1).A() * x_0; 
                r_norm = norm2(r);
                rel_norm = r_norm/r0_norm; 

                // print iteration status on every iteration 
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm}); 

                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1); 
                it++; 
            
            }


#ifdef BENCHMARKING_mode
                this->verbose(true); 
                it = 0; 
                converged = false; 
                r_norm = 1e11; 
                Scalar dr_norm, rho, rho_L2, e_norm, energy, corr, grid_complexity = 0.0, operator_complexity = 0.0, nnz_o = 0.0, nnz_g = 0.0; 
                Vector x_exact = x_0, e, Dr, x_k1, x_k01; 
                x_0 = 0. * x_exact; 
                D_inv = diag(level(l-1).A()); 
                D_inv = 1./D_inv; 
                
                Matrix A            =  level(l-1).A(); 
                Matrix P            =  transfers(l-2).I(); 
                
                grid_complexity     = get_global_nnz(P); 
                operator_complexity = get_global_nnz(A); 
                
                nnz_o                = operator_complexity; 
                nnz_g                = grid_complexity; 

                for(SizeType i = (l-2); i >= 0 ; i --)
                {                    
                    A =  level(i).A(); 
                    operator_complexity += get_global_nnz(A); 
                    
                    if(i <= l-3)
                    {
                        P =  transfers(i).I(); 
                        grid_complexity += get_global_nnz(P); 
                    }
                }

                operator_complexity /= nnz_o; 
                grid_complexity /= nnz_g; 

                this->init_solver("Multigrid - benchmarking \n C_op:   " + std::to_string(operator_complexity) + "\n C_gr:   " + std::to_string(grid_complexity) + + "\n levels:   " + std::to_string(l) , {" it. ", "|| r_N ||", "r_norm", "D^{-1} r", "|| e ||", "||e||_A - e^* x", "rho_L2", "rho_energy", "||corr||"}); 
                PrintInfo::print_iter_status(it, {r_norm});
                it++; 

                while(!converged)
                {        
                    if(it == 1)
                        x_k01 = x_0; 
                    x_k1 = x_0; 

                    multiplicative_cycle(rhs, l, x_0); 

                    r = rhs - level(l-1).A() * x_0; 
                    r_norm = norm2(r);
                    rel_norm = r_norm/r0_norm; 
                    dr_norm = dot(D_inv, r);
                    e = x_exact - x_0; 
                    e_norm  = norm2(e); 
                    energy = 0.5 * dot(e, level(l-1).A() *e)  + dot(e,x_0); 

                    if(it>1)
                    {
                        Vector next =  x_0 - x_k1; 
                        Vector prev = x_k01 - x_0;

                        rho_L2      =  Scalar(norm2(next))/ Scalar(norm2(prev)); 
                        rho         =  dot(next, level(l-1).A() *next)/ dot(prev, level(l-1).A() *prev); 
                        x_k01 = x_k1; 
                    }

                    corr = norm2(x_k1 - x_0); 
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm, dr_norm, e_norm, energy, rho_L2, rho, corr}); 

                    // check convergence and print interation info
                    converged = this->check_convergence(it, r_norm, rel_norm, 1); 
                    it++; 
                
                }
#endif

            return true; 
        }

        /**
         * @brief      The solve function for linear MG. 
         *             Needed in case, you want to use MG as LS in Nonlinear solvers.  
         *
         * @param[in]  A     The A.
         * @param[in]  b     The rhs.
         * @param      x0    The initial guess/solution.
         *
         * @return     
         */
        virtual bool solve(const Matrix &A, const Vector &b, Vector &x0) override
        {   
            update(make_ref(A));
            return solve(b, x0); 
        }

        inline Level &level(const SizeType &l)
        {
            return this->_levels[l]; 
        }

/*=======================================================================================================================================        =
=========================================================================================================================================*/
    private:
      

        inline Transfer &transfers(const SizeType & l)
        {
            return this->_transfers[l]; 
        }

        /**
         * @brief      Function implements multiplicative multigrid cycle. 
         *
         * @param[in]  rhs   The rhs.
         * @param[in]  l     The level.
         * @param      x_0   The current iterate. 
         *
         */
        virtual bool multiplicative_cycle(const Vector &rhs, const SizeType &l, Vector &x_0)
        {            
            assert(memory.valid(this->num_levels()) && l <= this->num_levels());
           
            Vector &r_h = memory.r_h[l-1];
            Vector &r_H = memory.r_H[l-1];
            Vector &c_H = memory.c_H[l-1];
            Vector &c_h = memory.c_h[l-1];

            if(local_size(x_0).get(0) != local_size(rhs).get(0)) {
               assert( local_size(x_0).get(0) != local_size(rhs).get(0) );
               std::cerr << "wrong local size for x_0 (if needed use redistribute_as(x_0, rhs)." << std::endl;
            }

            //SHOULD NOT BE NECEASSARY HERE
            // if(this->num_levels() > 2) {
            //     x_0 = redistribute_as(x_0, rhs); 
            // }

            // presmoothing 
            smoothing(level(l-1).A(), rhs, x_0, this->pre_smoothing_steps()); 

            // residual transfer 
            r_h = rhs - level(l-1).A() * x_0; 
            transfers(l-2).restrict(r_h, r_H); 

            // prepare correction 
            if(empty(c_H) || size(c_H).get(0) != size(r_H).get(0)) {
                c_H = local_zeros(local_size(r_H).get(0));        
            } else {
                c_H.set(0.);
            }

            if(l == 2) {
                // coarse solve 
                coarse_solve(level(l-2).A(), r_H, c_H);         

            } else {
                // recursive call into mg
                for(SizeType k = 0; k < this->mg_type(); k++) {   
                    SizeType l_new = l - 1; 
                    multiplicative_cycle(r_H, l_new, c_H); 
                }
            }

            // correction transfer
            transfers(l-2).interpolate(c_H, c_h); 

            if(use_line_search_) {
                const Scalar alpha = dot(c_h, r_h)/dot(level(l-1).A() * c_h, c_h);

                if(alpha <= 0.) {
                    // std::cerr << l << " zero grid correction" << std::endl;
                    x_0 += c_h;
                } else {
                    // std::cout << l << " : " << alpha << std::endl;
                    x_0 += alpha * c_h;
                }

            } else {
                x_0 += c_h; 
            }

            // postsmoothing 
            smoothing(level(l-1).A(), rhs, x_0, this->post_smoothing_steps()); 
            return true; 
        }

        /**
         * @brief      Function implements full multigrid cycle. 
         *              TODO:: fix
         *              - can be used jsut with homegenous BC - due to restriction of RHS  
         *              
         *              
         * @param[in]  rhs   The rhs.
         * @param[in]  l     The level.
         * @param      x_0   The current iterate. 
         *
         */
        virtual bool full_cycle(const Vector &rhs, const SizeType &l, Vector &x_0)
        {
            Vector rhs_h = rhs; 
            std::vector<Vector> rhss; 
            rhss.push_back(rhs);

            for(SizeType i = l-2; i >=0; i--)
            {
                transfers(i).restrict(rhs_h, rhs_h); 
                rhss.push_back(std::move(rhs_h));
            }

            coarse_solve(level(0).A(), rhss[l-1], x_0);    
            transfers(0).interpolate(x_0, x_0); 

            for(SizeType i = 1; i <l-1; i++)
            {
                for(SizeType j = 0; j < this->v_cycle_repetition(); j++) {
                    multiplicative_cycle(rhss[i], i+1, x_0);  
                }

                transfers(i).interpolate(x_0, x_0); 
            }

            for(SizeType i = 0; i < this->v_cycle_repetition(); i++) {
                multiplicative_cycle(rhss[0], l, x_0);
            }

            return true; 
        }

        /**
         * @brief      The function invokes smoothing. 
         *
         * @param[in]  A     The stifness matrix. 
         * @param[in]  rhs   The right hand side.
         * @param      x     The current iterate. 
         * @param[in]  nu    The number of smoothing steps. 
         *
         * @return   
         */
        bool smoothing(const Matrix &A, const Vector &rhs, Vector &x, const SizeType &nu = 1)
        {
            _smoother->sweeps(nu); 
            _smoother->smooth(A, rhs, x);

            return true; 
        }

        /**
         * @brief      The functions invokes coarse solve. 
         *
         * @param[in]  A    The A.
         * @param[in]  rhs  The right hand side.
         * @param      x    The current iterate.
         *
         * @return     
         */
        bool coarse_solve(const Matrix &A, const Vector &rhs, Vector &x)
        {
            _direct_solver->solve(A, rhs, x);
            return true; 
        }


    public:
        /**
         * @brief      Function changes direct solver needed for coarse grid solve. 
         *
         * @param[in]  linear_solver  The linear solver.
         *
         * @return     
         */
        bool change_direct_solver(const std::shared_ptr<Solver> &linear_solver = std::shared_ptr<Solver>())
        {
            _direct_solver = linear_solver; 
            _direct_solver->set_parameters(_parameters); 
            return true; 
        }

        /**
         * @brief      Function changes soother in MG. 
         *
         * @param[in]  smoother  The smoother.
         *
         * @return    
         */
        bool change_smoother(const std::shared_ptr<Smoother> &smoother = std::shared_ptr<Smoother>())
        {
            _smoother = smoother; 
            _smoother->set_parameters(_parameters); 
            return true; 
        }

        void set_perform_galerkin_assembly(const bool val)
        {
            perform_galerkin_assembly_ = val;   
        }

        void set_use_line_search(const bool val)
        {
            use_line_search_ = val;
        }

    protected:   
        std::shared_ptr<Smoother>           _smoother;
        std::shared_ptr<Solver>             _direct_solver;

    private:
        Parameters                          _parameters;
        bool perform_galerkin_assembly_;
        bool use_line_search_;

    };

}

#endif //UTOPIA_MULTIGRID_HPP


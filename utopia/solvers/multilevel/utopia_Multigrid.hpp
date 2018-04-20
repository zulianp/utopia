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
#include "utopia_LinearMultiLevel.hpp"
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
    class Multigrid : public LinearMultiLevel<Matrix, Vector>,
    public IterativeSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        
        typedef utopia::LinearSolver<Matrix, Vector>        Solver;
        typedef utopia::IterativeSolver<Matrix, Vector>     IterativeSolver;
        typedef utopia::Smoother<Matrix, Vector>            Smoother;
        typedef utopia::Level<Matrix, Vector>               Level;
        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef std::shared_ptr<Smoother> SmootherPtr;
        
        typedef struct
        {
            std::vector<Vector> r, c, c_I, r_R;
            
            void init(const int n_levels)
            {
                r.resize(n_levels);
                c.resize(n_levels);
                c_I.resize(n_levels);
                r_R.resize(n_levels);
            }
            
            inline bool valid(const std::size_t n_levels) const {
                return r.size() == n_levels;
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
         * @param[in]  coarse_solver  The direct solver for coarse level.
         */
        Multigrid(const std::shared_ptr<Smoother> &smoother,
                  const std::shared_ptr<Solver>   &coarse_solver,
                  const Parameters params = Parameters())
        : smoother_cloneable_(smoother),
          coarse_solver_(coarse_solver),
          perform_galerkin_assembly_(true),
          use_line_search_(false),
          block_size_(1)
        {
            set_parameters(params);
        }
        
        virtual ~Multigrid(){}
        
        void set_parameters(const Parameters params) override
        {
            IterativeSolver::set_parameters(params);
            LinearMultiLevel<Matrix, Vector>::set_parameters(params);
            smoother_cloneable_->set_parameters(params);
            
        }
        
        /*! @brief if overriden the subclass has to also call this one first
         */
        virtual void update(const std::shared_ptr<const Matrix> &op) override
        {
            IterativeSolver::update(op);
            
            if(perform_galerkin_assembly_) {
                this->galerkin_assembly(op);
            }
            
            smoothers_.resize(this->n_levels());
            smoothers_[0] = nullptr;
            
            for(std::size_t l = 1; l != smoothers_.size(); ++l) {
                smoothers_[l] = std::shared_ptr<Smoother>(smoother_cloneable_->clone());
                assert(smoothers_[l]);
                smoothers_[l]->update(level(l).A_ptr());
            }
            
            coarse_solver_->update(level(0).A_ptr());
        }
        
        /**
         * @brief      The solve function for multigrid method.
         *
         * @param[in]  rhs   The right hand side.
         * @param      x   The initial guess.
         *
         */
        virtual bool apply(const Vector &rhs, Vector &x) override
        {
            Scalar r_norm, r0_norm, rel_norm;
            SizeType it = 0;
            bool converged = false;
            SizeType L = this->n_levels();
            
            memory.init(L);
            SizeType l = L - 1;
            
            memory.r[l] = rhs - level(l).A() * x;
            
            r_norm = norm2(memory.r[l]);
            r0_norm = r_norm;
            
            this->init_solver("Multigrid", {" it. ", "|| r_N ||", "r_norm" });
            
            if(this->verbose())
                PrintInfo::print_iter_status(it, {r_norm, 1});
            it++;
            
            if(this->cycle_type() == FULL_CYCLE) {
                this->max_it(1);
            }
            
            while(!converged)
            {
                if(this->cycle_type() == MULTIPLICATIVE_CYCLE) {
                    multiplicative_cycle(l);
                } else if(this->cycle_type() == FULL_CYCLE) {
                    full_cycle(l);
                } else {
                    std::cout<<"ERROR::UTOPIA_MG<< unknown MG type... \n";
                }
                
#ifdef CHECK_NUM_PRECISION_mode
                if(has_nan_or_inf(x))
                {
                    x = local_zeros(local_size(x));
                    return true;
                }
#else
                // assert(!has_nan_or_inf(x));
#endif    
                
                x += memory.c[l];
                memory.r[l] = rhs - level(l).A() * x;
                r_norm = norm2(memory.r[l]);
                rel_norm = r_norm/r0_norm;
                
                // print iteration status on every iteration
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm});
                
                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1);
                it++;
                
            }
            
            
            return true;
        }
        
        inline Level &level(const SizeType &l)
        {
            return this->levels_[l];
        }
        
        /*=======================================================================================================================================        =
         =========================================================================================================================================*/
    private:
        
        
        /**
         * @brief      Function implements multiplicative multigrid cycle.
         *
         * @param[in]  l     The level.
         *
         */
        virtual bool multiplicative_cycle(const SizeType &l)
        {
            assert(memory.valid(this->n_levels()) && l < this->n_levels());
            
            Vector &r   = memory.r[l];
            Vector &c   = memory.c[l];
            Vector &c_I = memory.c_I[l];
            Vector &r_R = memory.r_R[l];
            
            if(empty(c)) {
                c = local_zeros(local_size(r).get(0));
            } else {
                c.set(0.);
            }

            assert(approxeq(double(norm2(c)), 0.));
            
            ////////////////////////////////////
            if(l == 0) {
                return coarse_solve(r, c);
            }
            ////////////////////////////////////
            
            // recursive call into mg
            for(SizeType k = 0; k < this->mg_type(); k++) {
                // presmoothing
                smoothing(l, r, c, this->pre_smoothing_steps());
                
                r_R = r - level(l).A() * c;
                
                // residual transfer
                this->transfer(l-1).restrict(r_R, memory.r[l-1]);
                assert(!empty(memory.r[l-1]));
                
                
                multiplicative_cycle(l-1);
                
                // correction transfer
                this->transfer(l-1).interpolate(memory.c[l-1], c_I);

#ifndef NDEBUG
                const Scalar err = norm2(r_R);
#endif
                if(use_line_search_) {
                    const Scalar alpha = dot(c_I, r_R)/dot(level(l).A() * c_I, c_I);
                    
                    if(alpha <= 0.) {
                        // std::cerr << l << " zero grid correction" << std::endl;
                        c += c_I;
                    } else {
                        // std::cout << l << " : " << alpha << std::endl;
                        c += alpha * c_I;
                    }
                    
                } else {
                    c += c_I;
                }

#ifndef NDEBUG
                const Scalar new_err = norm2(r - level(l).A() * c);
                // assert(new_err < err)
                if(new_err > err) {
                  std::cerr << "[Error] Multigrid::multiplicative: " << new_err << "<" << err << std::endl;
                }
#endif
                
                // postsmoothing
                smoothing(l, r, c, this->post_smoothing_steps());
            }
            
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
         * @param      x   The current iterate.
         *
         */
        virtual bool full_cycle(const SizeType &l)
        {
            Vector &r   = memory.r[l];
            Vector &c   = memory.c[l];
            Vector &c_I = memory.c_I[l];
            Vector &r_R = memory.r_R[l];
          
            std::vector<Vector> rhss(this->n_levels());
            rhss[l] = r;
            
            for(SizeType i = l - 1; i >= 0; i--)
            {
                this->transfer(i).restrict(rhss[i + 1], rhss[i]);
            }
            
            coarse_solve(rhss[0], memory.c[0]);
            this->transfer(0).interpolate(memory.c[0], memory.c[1]);
            
            for(SizeType i = 1; i < l; i++)
            {
                for(SizeType j = 0; j < this->v_cycle_repetition(); j++) {
                    memory.r[i] = rhss[i];
                    multiplicative_cycle(i);
                }
                
                this->transfer(i).interpolate(memory.c[i], memory.c[i+1]);
            }
            
            for(SizeType i = 0; i < this->v_cycle_repetition(); i++) {
                // memory.r[l] = rhss[l];
                multiplicative_cycle(l);
            }
            
            return true;
        }
        
        /**
         * @brief      The function invokes smoothing.
         *
         * @param[in]  level The level we are at.
         * @param[in]  rhs   The right hand side.
         * @param      x     The current iterate.
         * @param[in]  nu    The number of smoothing steps.
         *
         * @return
         */
        inline bool smoothing(const SizeType l, const Vector &rhs, Vector &x, const SizeType &nu = 1)
        {
            smoothers_[l]->sweeps(nu);
            smoothers_[l]->smooth(rhs, x);
            return true;
        }
        
        /**
         * @brief      The functions invokes coarse solve.
         *
         * @param[in]  rhs  The right hand side.
         * @param      x    The current iterate.
         *
         * @return
         */
        bool coarse_solve(const Vector &rhs, Vector &x)
        {
            if(!coarse_solver_->apply(rhs, x)) return false;
            assert(approxeq(level(0).A() * x, rhs, 1e-6));
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
        bool change_coarse_solver(const std::shared_ptr<Solver> &coarse_solver)
        {
            coarse_solver_ = coarse_solver;
            return true;
        }
        
        /**
         * @brief      Function changes soother in MG.
         *
         * @param[in]  smoother  The smoother.
         *
         * @return
         */
        bool change_smoother(const std::shared_ptr<Smoother> &smoother)
        {
            smoother_cloneable_ = smoother;
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
        
        inline void block_size(const int size)
        {
            block_size_ = size;
        }
        
    protected:
        std::shared_ptr<Smoother> smoother_cloneable_;
        std::shared_ptr<Solver>   coarse_solver_;
        std::vector<SmootherPtr>  smoothers_;
        
    private:
        bool perform_galerkin_assembly_;
        bool use_line_search_;
        int block_size_;
        
    };
    
}

#endif //UTOPIA_MULTIGRID_HPP


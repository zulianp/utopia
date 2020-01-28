#ifndef UTOPIA_MULTIGRID_QR_HPP
#define UTOPIA_MULTIGRID_QR_HPP
#include "utopia_Smoother.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Core.hpp"
#include "utopia_Utils.hpp"
#include "utopia_PrintInfo.hpp"
#include "utopia_ConvergenceReason.hpp"
#include "utopia_Level.hpp"
#include "utopia_LinearMultiLevel.hpp"
#include "utopia_ProjectedGaussSeidelQR.hpp"
#include "utopia_Recorder.hpp"
#include "utopia_MatrixTruncatedTransfer.hpp"

#include <ctime>
#include <cassert>

namespace utopia
{
    /**
     * @brief      The class for Linear Multigrid solver.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class MultigridQR final:  public LinearMultiLevel<Matrix, Vector>,
                              public IterativeSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef utopia::LinearSolver<Matrix, Vector>        Solver;
        typedef utopia::IterativeSolver<Matrix, Vector>     IterativeSolver;
        typedef utopia::Smoother<Matrix, Vector>            Smoother;
        typedef utopia::Level<Matrix, Vector>               Level;
        typedef utopia::Transfer<Matrix, Vector>            Transfer;
        typedef std::shared_ptr<Smoother>                   SmootherPtr;

        typedef struct
        {
            std::vector<Vector> r, c, x, rhs;

            void init(const int n_levels)
            {
                r.resize(n_levels);
                c.resize(n_levels);
                x.resize(n_levels);
                rhs.resize(n_levels);
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
        MultigridQR(const std::shared_ptr<Smoother> &smoother,
                  const std::shared_ptr<Solver>   &coarse_solver, 
                  const SizeType & num_levels)
        : smoother_cloneable_(smoother),
          coarse_solver_(coarse_solver),
          perform_galerkin_assembly_(true),
          use_line_search_(false),
          block_size_(1)
        {
            this->must_generate_masks(true);
            this->num_levels_ = num_levels; 

            smoothers_.resize(this->n_levels());   
        }

        ~MultigridQR(){}

        void read(Input &in) override
        {
          LinearMultiLevel<Matrix, Vector>::read(in);
          IterativeSolver::read(in);

          in.get("perform_galerkin_assembly", perform_galerkin_assembly_);
          in.get("use_line_search", use_line_search_);
          in.get("block_size", block_size_);

          if(smoother_cloneable_) {
              in.get("smoother", *smoother_cloneable_);
          }
          if(coarse_solver_) {
              in.get("coarse_solver", *coarse_solver_);
          }

        }

        void print_usage(std::ostream &os) const override
        {
          LinearMultiLevel<Matrix, Vector>::print_usage(os);
          IterativeSolver::print_usage(os);

          this->print_param_usage(os, "perform_galerkin_assembly", "bool", "Flag turning on/off galerkin assembly.", "true");
          this->print_param_usage(os, "use_line_search", "bool", "Flag turning on/off line-search after coarse grid correction.", "false");
          this->print_param_usage(os, "block_size", "int", "Block size for systems.", "1");
          this->print_param_usage(os, "smoother", "Smoother", "Input parameters for all smoothers.", "-");
          this->print_param_usage(os, "coarse_solver", "LinearSolver", "Input parameters for coarse solver.", "-");
        }

        /*! @brief if overriden the subclass has to also call this one first
         */
        void update(const std::shared_ptr<const Matrix> &op) override
        {
            IterativeSolver::update(op);

            if(perform_galerkin_assembly_) {
                this->galerkin_assembly(op);
            }

            update();
        }

        /*! @brief if no galerkin assembly is used but instead set_linear_operators is used.
                   One can call this update instead of the other one.
         */
        void update()
        {
            smoothers_.resize(this->n_levels());
            smoothers_[0] = nullptr;

            for(std::size_t l = 1; l != smoothers_.size(); ++l) {
              if(smoothers_[l]==nullptr){
                  smoothers_[l] = std::shared_ptr<Smoother>(smoother_cloneable_->clone());
                  assert(smoothers_[l]);
                  std::cout<<"UTOPIA_QR::smoothers got reseted.... \n"; 
                }
                smoothers_[l]->update(level(l).A_ptr());
            }

            coarse_solver_->update(level(0).A_ptr());
        }

        // rememebr to put level -1, as c++ indexing goes from 0
        void set_smoother(const std::shared_ptr<Smoother> &smoother, const SizeType & level)
        {
          if(level > this->n_levels() || level < 1){
            utopia_error("MultigridQR:: set_smoother, invalid level.");
          }

          if(smoothers_.size()-1 < level){
            utopia_error("MultigridQR:: set_smoother, smoothers array was not allocated.");
          }
          smoothers_[level] = smoother; 
        }

        /**
         * @brief      The solve function for multigrid method.
         *
         * @param[in]  rhs   The right hand side.
         * @param      x   The initial guess.
         *
         */
        bool apply(const Vector &rhs, Vector &x_fine) override
        {

            // UTOPIA_RECORD_SCOPE_BEGIN("apply");
            Scalar r_norm, r0_norm, rel_norm;
            SizeType it = 0;
            bool converged = false;
            bool ok = true; UTOPIA_UNUSED(ok);
            SizeType L = this->n_levels();

            memory.init(L);
            SizeType l = L - 1;

            // UTOPIA_RECORD_VALUE("rhs", rhs);
            // UTOPIA_RECORD_VALUE("x", x);
            memory.x[l] = x_fine;
            // memory.r[l] = rhs - level(l).A() * memory.x[l];
            memory.rhs[l] = rhs; 

            // UTOPIA_RECORD_VALUE("rhs - level(l).A() * x", memory.r[l]);

            r_norm = norm2(rhs - level(l).A() * memory.x[l]);
            r0_norm = r_norm;

            Vector x_old = memory.x[l]; 

            std::string mg_header_message = "Multigrid: " + std::to_string(L) +  " levels";
            this->init_solver(mg_header_message, {" it. ", "|| r_N ||", "||x_old - x_new||" });

            if(this->verbose())
                PrintInfo::print_iter_status(it, {r_norm, 1});
            it++;

            if(this->cycle_type() == FULL_CYCLE) {
                this->max_it(1);
            }

            while(!converged)
            {
                if(this->cycle_type() == MULTIPLICATIVE_CYCLE) {
                    ok = standard_cycle(l); assert(ok);
                } 
                else {
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

                // UTOPIA_RECORD_VALUE("memory.c[l]", memory.c[l]);

                //x += memory.c[l];

                rel_norm = norm2(x_old - memory.x[l]); 
                x_old = memory.x[l]; 

                // UTOPIA_RECORD_VALUE("x += memory.c[l]", x);

                // memory.r[l] = rhs - level(l).A() * memory.x[l];
                r_norm = norm2(rhs - level(l).A() * memory.x[l]);
                // rel_norm = r_norm/r0_norm;

                // print iteration status on every iteration
                if(this->verbose())
                    PrintInfo::print_iter_status(it, {r_norm, rel_norm});

                // check convergence and print interation info
                converged = this->check_convergence(it, r_norm, rel_norm, 1);
                it++;

            }

            this->print_statistics(it);

            // UTOPIA_RECORD_SCOPE_END("apply");
            return true;
        }

        inline Level &level(const SizeType &l)
        {
            return this->levels_[l];
        }

        /*=======================================================================================================================================        =
         =========================================================================================================================================*/
    private:
        bool standard_cycle(const SizeType &l)
        {
            // UTOPIA_RECORD_SCOPE_BEGIN("standard_cycle(" + std::to_string(l) + ")");
            assert(memory.valid(this->n_levels()) && l < this->n_levels());

            // Vector &r   = memory.r[l];
            // Vector &c   = memory.c[l];
            // Vector &x   = memory.x[l];
            // Vector &rhs = memory.rhs[l];

            // if(empty(c) || size(c) != size(r)) {
            //     c = local_zeros(local_size(r).get(0));
            // } else {
            //     c.set(0.);
            // }

            // assert(approxeq(double(norm2(c)), 0.));

            ////////////////////////////////////
            if(l == 0) 
            {
                coarse_solve(memory.r[l], memory.c[l]);

                // coarse_solve(r, c);
            //     // if(coarse_solve(r, c)) 
            //     // {
            //     //   assert(approxeq(level(0).A() * c, r, 1e-6));
            //     //   // UTOPIA_RECORD_VALUE("coarse_solve(r, c)", c);
            //     //   // UTOPIA_RECORD_SCOPE_END("standard_cycle(" + std::to_string(l) + ")");
            //       return true;
            //     // } 
            //     // else {
            //     //   assert(false);
                   return false;
            //     // }
            }
            ////////////////////////////////////

            // recursive call into mg
            for(SizeType k = 0; k < this->mg_type(); k++) {
                // presmoothing

                if(l == this->n_levels()-1)
                {
                    // do fine level smoothing
                    smoothing_fine(l, memory.rhs[l], memory.x[l], this->pre_smoothing_steps(), true);
                    memory.r[l] = memory.rhs[l] - level(l).A() * memory.x[l];
                }
                else if(l > 0)
                {
                    smoothing(l, memory.r[l], memory.c[l], this->pre_smoothing_steps());
                }

                // UTOPIA_RECORD_VALUE("smoothing(l, r, c, this->pre_smoothing_steps());", c);


                // memory.rhs[l] = r - level(l).A() * x;
                // residual transfer
                this->transfer(l-1).restrict(memory.r[l], memory.r[l-1]);
                // memory.rhs[l-1] = memory.r[l-1];

                // UTOPIA_RECORD_VALUE("this->transfer(l-1).restrict(r_R, memory.r[l-1]);", memory.r[l-1]);

                //NEW
                // if(this->must_generate_masks()) {
                //   this->apply_mask(l-1, memory.rhs[l-1]);
                //   // UTOPIA_RECORD_VALUE("this->apply_mask(l-1, memory.r[l-1]);", memory.r[l-1]);
                // }

                // assert(!empty(memory.rhs[l-1]));

                memory.c[l-1] = 0*memory.r[l-1]; 
                standard_cycle(l-1);
                
                // correction transfer
                this->transfer(l-1).interpolate(memory.c[l-1], memory.c[l]);

                memory.x[l] += memory.c[l];

                // postsmoothing
                if( l == this-> n_levels()-1)
                {
                  smoothing_fine(l, memory.rhs[l], memory.x[l], this->post_smoothing_steps(), false);
                  memory.r[l] = memory.rhs[l] - level(l).A() * memory.x[l];

                }
                else if(l > 0)
                {
                  smoothing(l, memory.r[l-1], memory.c[l-1], this->post_smoothing_steps());
                  //std::cerr<< "--------- not fixed yet.... \n"; 
                }
                // UTOPIA_RECORD_VALUE("smoothing(l, r, c, this->post_smoothing_steps());", c);

            }

            // UTOPIA_RECORD_SCOPE_END("standard_cycle(" + std::to_string(l) + ")");
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
            // GaussSeidel<Matrix, Vector>* GS_smoother =  dynamic_cast<GaussSeidel<Matrix, Vector>* > (smoothers_[l].get()); 

            std::cout<<"----- regular smoothing --- "<< std::endl;

            smoothers_[l]->sweeps(nu);
            smoothers_[l]->smooth(rhs, x);

            return true;
        }

        
        inline bool smoothing_fine(const SizeType l, const Vector &rhs, Vector &x, const SizeType &nu = 1, const bool & pre_sm = false)
        {
            assert(l != this->n_levels()-2); 
            
            if(ProjectedGaussSeidelQR<Matrix, Vector>* GS_smoother =  dynamic_cast<ProjectedGaussSeidelQR<Matrix, Vector>* > (smoothers_[l].get())){
              GS_smoother->verbose(true);
              GS_smoother->set_R(R_);            
              GS_smoother->sweeps(nu);
              GS_smoother->set_box_constraints(make_box_constaints(make_ref(lb_),  make_ref(ub_)));          
              GS_smoother->smooth(rhs, x);
            }
            else{
              utopia_error("MG_QR: requires ProjectedGaussSeidelQR to be the fine level smoother "); 
            }
            if(pre_sm){
             
             if(MatrixTruncatedTransfer<Matrix, Vector>* trunc_transfer =  dynamic_cast<MatrixTruncatedTransfer<Matrix, Vector>* > (this->transfers_[l-1].get())){ 
                ProjectedGaussSeidelQR<Matrix, Vector>* GS_smoother =  dynamic_cast<ProjectedGaussSeidelQR<Matrix, Vector>* > (smoothers_[l].get());
                const Vector & inactive_set = GS_smoother->get_inactive_set();
                trunc_transfer->truncate_interpolation(inactive_set); 
                this->galerkin_assembly(this->get_operator());
                // update();
                for(std::size_t l = 1; l != smoothers_.size()-1; ++l) {
                    smoothers_[l]->update(level(l).A_ptr());
                }
                coarse_solver_->update(level(0).A_ptr());              
              }
            else{
              utopia_error("MG_QR: requires MatrixTruncatedTransfer for the finest level.");
            }
          }
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

        MultigridQR * clone() const override
        {
           return new MultigridQR(
            std::shared_ptr<Smoother>(smoother_cloneable_->clone()),
            std::shared_ptr<Solver>(coarse_solver_->clone()), 
            this->n_levels()
            );
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

        void use_line_search(const bool val)
        {
            use_line_search_ = val;
        }

        bool use_line_search() const
        {
            return use_line_search_;
        }

        inline void block_size(const SizeType size)
        {
            block_size_ = size;
        }

        inline SizeType block_size() const
        {
          return block_size_;
        }

        void set_R(const Matrix & R)
        {
          R_ = R; 
        }

        void set_lower_bound(const Vector & lb)
        {
          lb_ = lb; 
        }        

        void set_upper_bound(const Vector & ub)
        {
          ub_ = ub; 
        }                
        
        // virtual bool set_transfer_operators(const std::vector<std::shared_ptr<Matrix>> &interpolation_operators)
        // {
        //     if(this->n_levels() <= 0){
        //         this->n_levels(interpolation_operators.size() + 1);
        //     }
        //     else if(this->n_levels() != static_cast<SizeType>(interpolation_operators.size()) + 1){
        //         utopia_error("utopia::MultilevelBase:: number of levels and transfer operators do not match ... \n");
        //     }

        //     this->transfers_.clear();
        //     for(auto I = interpolation_operators.begin(); I != interpolation_operators.end() ; ++I )
        //         this->transfers_.push_back(std::make_shared<MatrixTransfer>(*I));

        //     return true;
        // }

    protected:
        std::shared_ptr<Smoother> smoother_cloneable_;
        std::shared_ptr<Solver>   coarse_solver_;
        std::vector<SmootherPtr>  smoothers_;


    private:
        bool perform_galerkin_assembly_;
        bool use_line_search_;
        SizeType block_size_;
        Matrix Q_, R_;
        Vector lb_, ub_; 
    };

}

#endif //UTOPIA_MULTIGRID_QR_HPP


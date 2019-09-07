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


#include "utopia_Recorder.hpp"

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
                  const std::shared_ptr<Solver>   &coarse_solver)
        : smoother_cloneable_(smoother),
          coarse_solver_(coarse_solver),
          perform_galerkin_assembly_(true),
          use_line_search_(false),
          block_size_(1)
        {
            this->must_generate_masks(true);
        }

        ~Multigrid(){}

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
        bool apply(const Vector &rhs, Vector &x) override
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

            memory.r[l] = rhs - level(l).A() * x;

            // UTOPIA_RECORD_VALUE("rhs - level(l).A() * x", memory.r[l]);

            r_norm = norm2(memory.r[l]);
            r0_norm = r_norm;

            std::string mg_header_message = "Multigrid: " + std::to_string(L) +  " levels";
            this->init_solver(mg_header_message, {" it. ", "|| r_N ||", "r_norm" });

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
                } else if(this->cycle_type() == FULL_CYCLE) {
                    ok = full_cycle(l); assert(ok);
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

                // UTOPIA_RECORD_VALUE("memory.c[l]", memory.c[l]);

                x += memory.c[l];

                // UTOPIA_RECORD_VALUE("x += memory.c[l]", x);

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


        /**
         * @brief      Function implements multigrid standard cycle.
         *
         * @param[in] l The level.
         *
         */
        bool standard_cycle(const SizeType &l)
        {
            // UTOPIA_RECORD_SCOPE_BEGIN("standard_cycle(" + std::to_string(l) + ")");
            assert(memory.valid(this->n_levels()) && l < this->n_levels());

            Vector &r   = memory.r[l];
            Vector &c   = memory.c[l];
            Vector &c_I = memory.c_I[l];
            Vector &r_R = memory.r_R[l];

            if(empty(c) || size(c) != size(r)) {
                c = local_zeros(local_size(r).get(0));
            } else {
                c.set(0.);
            }

            assert(approxeq(double(norm2(c)), 0.));

            ////////////////////////////////////
            if(l == 0) {
              // UTOPIA_RECORD_VALUE("c", c);
              // UTOPIA_RECORD_VALUE("r", r);
                if(coarse_solve(r, c)) {
                  assert(approxeq(level(0).A() * c, r, 1e-6));
                  // UTOPIA_RECORD_VALUE("coarse_solve(r, c)", c);
                  // UTOPIA_RECORD_SCOPE_END("standard_cycle(" + std::to_string(l) + ")");
                  return true;
                } else {
                  assert(false);
                  return false;
                }
            }
            ////////////////////////////////////

            // recursive call into mg
            for(SizeType k = 0; k < this->mg_type(); k++) {
                // presmoothing
                smoothing(l, r, c, this->pre_smoothing_steps());
                // UTOPIA_RECORD_VALUE("smoothing(l, r, c, this->pre_smoothing_steps());", c);


                r_R = r - level(l).A() * c;

                // UTOPIA_RECORD_VALUE("r_R", r_R);

                // residual transfer
                this->transfer(l-1).restrict(r_R, memory.r[l-1]);

                // UTOPIA_RECORD_VALUE("this->transfer(l-1).restrict(r_R, memory.r[l-1]);", memory.r[l-1]);

                //NEW
                if(this->must_generate_masks()) {
                  this->apply_mask(l-1, memory.r[l-1]);
                  // UTOPIA_RECORD_VALUE("this->apply_mask(l-1, memory.r[l-1]);", memory.r[l-1]);
                }

                assert(!empty(memory.r[l-1]));




                standard_cycle(l-1);

                // correction transfer
                this->transfer(l-1).interpolate(memory.c[l-1], c_I);
                // UTOPIA_RECORD_VALUE("this->transfer(l-1).interpolate(memory.c[l-1], c_I);", c_I);

#ifndef NDEBUG
                const Scalar err = norm2(r_R);
#endif
                if(use_line_search_) {
                    const Scalar alpha = dot(c_I, r_R)/dot(level(l).A() * c_I, c_I);

                    if(alpha <= 0.) {
                        std::cerr << l << " alpha: " << alpha << std::endl;
                        c += c_I;
                    } else {
                        // std::cout << l << " : " << alpha << std::endl;
                        c += alpha * c_I;
                    }

                } else {
                    c += c_I;
                }

                // UTOPIA_RECORD_VALUE("c", c);

                // postsmoothing
                smoothing(l, r, c, this->post_smoothing_steps());

                // UTOPIA_RECORD_VALUE("smoothing(l, r, c, this->post_smoothing_steps());", c);

#ifndef NDEBUG
                const Scalar new_err = norm2(r - level(l).A() * c);
                // assert(new_err < err)
                if(new_err > err) {
                  m_utopia_error_once("[Error] Multigrid::standard_cycle (" + std::to_string(l) + "): coarse grid correction raises error " + std::to_string(new_err) + "<" + std::to_string(err));
                }
#endif
            }

            // UTOPIA_RECORD_SCOPE_END("standard_cycle(" + std::to_string(l) + ")");
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
        bool full_cycle(const SizeType &l)
        {
            Vector &r   = memory.r[l];
            // Vector &c   = memory.c[l];
            // Vector &c_I = memory.c_I[l];
            // Vector &r_R = memory.r_R[l];

            std::vector<Vector> rhss(this->n_levels());
            rhss[l] = r;

            for(SizeType i = l - 1; i >= 0; i--)
            {
                this->transfer(i).restrict(rhss[i + 1], rhss[i]);
                //NEW
                this->apply_mask(i, rhss[i]);
            }

            coarse_solve(rhss[0], memory.c[0]);
            this->transfer(0).interpolate(memory.c[0], memory.c[1]);

            for(SizeType i = 1; i < l; i++)
            {
                for(SizeType j = 0; j < this->v_cycle_repetition(); j++) {
                    memory.r[i] = rhss[i];
                    standard_cycle(i);
                }

                this->transfer(i).interpolate(memory.c[i], memory.c[i+1]);
            }

            for(SizeType i = 0; i < this->v_cycle_repetition(); i++) {
                // memory.r[l] = rhss[l];
                standard_cycle(l);
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

        Multigrid * clone() const override
        {
           return new Multigrid(
            std::shared_ptr<Smoother>(smoother_cloneable_->clone()),
            std::shared_ptr<Solver>(coarse_solver_->clone())
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


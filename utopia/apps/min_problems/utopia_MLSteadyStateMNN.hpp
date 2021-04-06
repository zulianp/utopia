#ifndef UTOPIA_ML_STEADY_STATE_MNN_HPP
#define UTOPIA_ML_STEADY_STATE_MNN_HPP

#include "utopia_BlockQPSolver.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_Input.hpp"
#include "utopia_MassMatrix.hpp"
#include "utopia_MonotoneMultigrid.hpp"
#include "utopia_Multilevel.hpp"
#include "utopia_ProjectedChebyshev3level.hpp"
#include "utopia_RedundantQPSolver.hpp"
#include "utopia_TestFunctions.hpp"
#include "utopia_make_unique.hpp"

#include <memory>

namespace utopia {

    // FIXME complete the overriding process
    template <class FunctionSpace, class ProblemType, class BCType, class ICType>
    class MLSteadyStateMNN final : public Configurable {
    public:
        using Matrix = typename FunctionSpace::Matrix;
        using Vector = typename FunctionSpace::Vector;
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;

        MLSteadyStateMNN(FunctionSpace &space_coarse)
            : init_(false),
              n_levels_(2),
              log_output_path_("monotone_mg_log_file.csv"),
              output_path_("mnmg_out"),
              save_output_(true) {
            spaces_.resize(2);
            spaces_[0] = make_ref(space_coarse);
        }

        void read(Input &in) override {
            in.get("log_output_path", log_output_path_);
            in.get("output_path", output_path_);
            in.get("n_levels", n_levels_);
            in.get("save_output", save_output_);

            init_ml_setup();

            fun_->read(in);
            BC_conditions_->read(in);
            IC_->read(in);

            in.get("solver", *multigrid_);
        }

        bool init_ml_setup() {
            if (n_levels_ < 2) {
                std::cerr << "n_levels must be at least 2" << std::endl;
                return false;
            }

            spaces_.resize(n_levels_);

            transfers_.resize(n_levels_ - 1);

            for (SizeType i = 1; i < n_levels_; ++i) {
                spaces_[i] = spaces_[i - 1]->uniform_refine();

                auto I = std::make_shared<Matrix>();
                spaces_[i - 1]->create_interpolation(*spaces_[i], *I);
                assert(!empty(*I));

                Matrix Iu;  // = *I;
                Iu.destroy();
                MatConvert(raw_type(*I), I->type_override(), MAT_INITIAL_MATRIX, &raw_type(Iu));

                transfers_[i - 1] =
                    std::make_shared<IPRTruncatedTransfer<Matrix, Vector>>(std::make_shared<Matrix>(Iu));
            }

            // initial conddition needs to be setup only on the finest level
            IC_ = std::make_shared<ICType>(*spaces_.back());
            fun_ = std::make_shared<ProblemType>(*spaces_.back());
            BC_conditions_ = std::make_shared<BCType>(*spaces_.back());

            //////////////////////////////////////////////// init solver
            //////////////////////////////////////////////////////

            if (!multigrid_) {
                // auto fine_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
                // auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();

                // auto fine_smoother = std::make_shared<ProjectedChebyshev3level<Matrix, Vector>>();
                // auto coarse_smoother = std::make_shared<Chebyshev3level<Matrix, Vector>>();

                auto fine_smoother = std::make_shared<ProjectedChebyshev3level<Matrix, Vector>>();
                auto coarse_smoother = std::make_shared<GaussSeidel<Matrix, Vector>>();

                // auto fine_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
                // auto coarse_smoother = std::make_shared<Chebyshev3level<Matrix, Vector>>();

                auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>();

                multigrid_ = std::make_shared<MonotoneMultigrid<Matrix, Vector>>(
                    fine_smoother, coarse_smoother, direct_solver, n_levels_);
            }

            multigrid_->set_transfer_operators(transfers_);
            multigrid_->verbose(true);

            init_ = true;

            return true;
        }

        FunctionSpace &fine_space() { return *spaces_.back(); }

        const FunctionSpace &fine_space() const { return *spaces_.back(); }

        std::shared_ptr<FunctionSpace> fine_space_ptr() { return spaces_.back(); }

        std::shared_ptr<const FunctionSpace> fine_space_ptr() const { return spaces_.back(); }

        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////

        void init_solution(Vector &solution) {
            spaces_.back()->create_vector(solution);
            rename("X", solution);

            IC_->init(solution);
            spaces_.back()->apply_constraints(solution);
        }

        void write_to_file(FunctionSpace & /*space*/, const Vector &solution) {
            if (save_output_) {
                spaces_.back()->write(this->output_path_ + ".vtr", solution);
                Utopia::instance().set("log_output_path", log_output_path_);
            }
        }

        void prepare_for_solve(Vector &solution) {
            BC_conditions_->emplace_BC();

            // update fine level solution  and constraint
            spaces_.back()->apply_constraints(solution);

            Vector &bc_flgs = fun_->get_eq_constrains_flg();
            Vector &bc_values = fun_->get_eq_constrains_values();

            spaces_.back()->apply_constraints(bc_values);
            spaces_.back()->build_constraints_markers(bc_flgs);
            fun_->init_constraint_indices();
        }

        void init(FunctionSpace &space, Vector &solution) {
            init_solution(solution);
            write_to_file(space, solution);
        }

        void run() {
            if (!init_) {
                init_ml_setup();
            }

            Vector solution, lower_bound, upper_bound;

            // init fine level spaces
            this->init(*spaces_[n_levels_ - 1], solution);
            prepare_for_solve(solution);

            lower_bound = solution;
            upper_bound = solution;

            upper_bound.set(0.25);
            lower_bound.set(-0.25);

            solution.set(0.0);
            Vector g = 0.0 * solution;
            Matrix H;

            fun_->gradient(solution, g);
            fun_->hessian(solution, H);

            // multigrid_->set_box_constraints(make_lower_bound_constraints(make_ref(lower_bound)));
            multigrid_->set_box_constraints(make_box_constaints(make_ref(lower_bound), make_ref(upper_bound)));
            multigrid_->update(make_ref(H));

            multigrid_->active_set().tol(1e-15);
            multigrid_->max_it(40);
            multigrid_->atol(1e-9);
            multigrid_->apply(g, solution);

            write_to_file(*spaces_[n_levels_ - 1], solution);
        }

    private:
        bool init_;
        SizeType n_levels_;

        std::vector<std::shared_ptr<FunctionSpace>> spaces_;
        std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> transfers_;

        std::shared_ptr<ExtendedFunction<Matrix, Vector>> fun_;
        std::shared_ptr<BCType> BC_conditions_;
        std::shared_ptr<ICType> IC_;

        std::string log_output_path_;
        std::string output_path_;

        std::shared_ptr<MonotoneMultigrid<Matrix, Vector>> multigrid_;

        bool save_output_;
    };

}  // namespace utopia

#endif  // UTOPIA_ML_STEADY_STATE_MNN_HPP
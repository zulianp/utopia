#ifndef UTOPIA_DM_RMTR_SETUP_STEADY_STATE_HPP
#define UTOPIA_DM_RMTR_SETUP_STEADY_STATE_HPP

#include "utopia_BlockQPSolver.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_Input.hpp"
#include "utopia_MassMatrix.hpp"
#include "utopia_Multilevel.hpp"
#include "utopia_RedundantQPSolver.hpp"
#include "utopia_make_unique.hpp"

#include <memory>

namespace utopia {

// FIXME complete the overriding process
template <class FunctionSpace, class ProblemType, class BCType>
class MLSteadyState final : public Configurable {
 public:
  using Matrix = typename FunctionSpace::Matrix;
  using Vector = typename FunctionSpace::Vector;
  using Scalar = typename FunctionSpace::Scalar;
  using SizeType = typename FunctionSpace::SizeType;

  MLSteadyState(FunctionSpace &space_coarse)
      : init_(false),
        n_levels_(2),
        n_coarse_sub_comm_(1),
        log_output_path_("rmtr_log_file.csv"),
        output_path_("rmtr_out.csv"),
        save_output_(true){
    spaces_.resize(2);
    spaces_[0] = make_ref(space_coarse);
  }

  void read(Input &in) override {
  
    in.get("log_output_path", log_output_path_);
    in.get("output_path", output_path_);
    in.get("n_coarse_sub_comm", n_coarse_sub_comm_);
    in.get("n_levels", n_levels_);
    in.get("save_output", save_output_);

    init_ml_setup();

    for (std::size_t l = 0; l < level_functions_.size(); l++) {
      level_functions_[l]->read(in);
      BC_conditions_[l]->read(in);
    }

    in.get("solver", *rmtr_);
  }

  bool init_ml_setup() {
    if (n_levels_ < 2) {
      std::cerr << "n_levels must be at least 2" << std::endl;
      return false;
    }

    spaces_.resize(n_levels_);

    level_functions_.resize(n_levels_);
    auto fun = std::make_shared<ProblemType>(*spaces_[0]);
    // fun->use_crack_set_irreversibiblity(false);
    level_functions_[0] = fun;

    BC_conditions_.resize(n_levels_);
    auto bc = std::make_shared<BCType>(*spaces_[0]);
    BC_conditions_[0] = bc;

    transfers_.resize(n_levels_ - 1);

    for (SizeType i = 1; i < n_levels_; ++i) {
      spaces_[i] = spaces_[i - 1]->uniform_refine();

      auto fun = std::make_shared<ProblemType>(*spaces_[i]);

      level_functions_[i] = fun;

      auto bc = std::make_shared<BCType>(*spaces_[i]);
      BC_conditions_[i] = bc;

      auto I = std::make_shared<Matrix>();
      spaces_[i - 1]->create_interpolation(*spaces_[i], *I);
      assert(!empty(*I));

      Matrix Iu;  // = *I;
      Iu.destroy();
      MatConvert(raw_type(*I), I->type_override(), MAT_INITIAL_MATRIX,
                 &raw_type(Iu));
      Matrix R = transpose(Iu);

      Reaction<FunctionSpace> mass_matrix_assembler_fine(*spaces_[i]);
      Matrix M_fine;
      mass_matrix_assembler_fine.mass_matrix(M_fine);

      Reaction<FunctionSpace> mass_matrix_assembler_coarse(*spaces_[i - 1]);
      Matrix M_coarse;
      mass_matrix_assembler_coarse.mass_matrix(M_coarse);

      Matrix inv_lumped_mass = diag(1. / sum(M_coarse, 1));
      Matrix P = inv_lumped_mass * R * M_fine;

      transfers_[i - 1] = std::make_shared<IPTransferNested<Matrix, Vector>>(
          std::make_shared<Matrix>(Iu), std::make_shared<Matrix>(P));
    }

    //////////////////////////////////////////////// init solver
    //////////////////////////////////////////////////////

    if (!rmtr_) {
      rmtr_ = std::make_shared<RMTR_inf<
          Matrix, Vector, TRBoundsGratton<Matrix, Vector>, SECOND_ORDER> >(
          n_levels_);
    }

    // std::shared_ptr<QPSolver<PetscMatrix, PetscVector>> tr_strategy_fine =
    auto tr_strategy_fine = std::make_shared<utopia::ProjectedGaussSeidel<Matrix, Vector>>();
    tr_strategy_fine->n_local_sweeps(1); 


    std::shared_ptr<QPSolver<Matrix, Vector>> tr_strategy_coarse =
        std::make_shared<SemismoothNewton<Matrix, Vector>>(
            std::make_shared<Factorization<Matrix, Vector>>());

    rmtr_->verbosity_level(utopia::VERBOSITY_LEVEL_VERY_VERBOSE);
    // rmtr_->verbosity_level(utopia::VERBOSITY_LEVEL_NORMAL);

    rmtr_->set_coarse_tr_strategy(tr_strategy_coarse);
    rmtr_->set_fine_tr_strategy(tr_strategy_fine);

    rmtr_->set_transfer_operators(transfers_);
    rmtr_->set_functions(level_functions_);
    rmtr_->verbose(true);
    rmtr_->atol(1e-9);


    init_ = true;

    return true;
  }

  FunctionSpace &fine_space() { return *spaces_.back(); }

  const FunctionSpace &fine_space() const { return *spaces_.back(); }

  std::shared_ptr<FunctionSpace> fine_space_ptr() { return spaces_.back(); }

  std::shared_ptr<const FunctionSpace> fine_space_ptr() const {
    return spaces_.back();
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////


  void write_to_file(FunctionSpace &space, const Vector &solution) {
    if (save_output_) {
      spaces_.back()->write(this->output_path_ + ".vtr", solution);
      Utopia::instance().set("log_output_path", log_output_path_);
    }
  }

  void prepare_for_solve(Vector &solution) {
    for (std::size_t l = 0; l < BC_conditions_.size(); l++) {
      BC_conditions_[l]->emplace_BC();
    }

    // update fine level solution  and constraint
    spaces_.back()->apply_constraints(solution);

    for (std::size_t l = 0; l < BC_conditions_.size(); l++) {
      Vector &bc_flgs = level_functions_[l]->get_eq_constrains_flg();
      Vector &bc_values = level_functions_[l]->get_eq_constrains_values();

      spaces_[l]->apply_constraints(bc_values);
      spaces_[l]->build_constraints_markers(bc_flgs);
      level_functions_[l]->init_constraint_indices();
    }
  }

  void init(FunctionSpace &space, Vector &solution) {
    write_to_file(space, solution);
  }

  void run() {
    if (!init_) {
      init_ml_setup();
    }

    // init fine level spaces
    Vector solution = level_functions_.back()->initial_guess(); 
    this->init(*spaces_[n_levels_ - 1], solution);



    prepare_for_solve(solution);


    
    rmtr_->delta0(1e9);
    rmtr_->solve(solution);
    auto sol_status = rmtr_->solution_status();
    sol_status.describe(std::cout); 
    // rmtr_->print_usage(std::cout); 



    // ////////////////////////////////////////////////////////////////////////////////////////////////////////
    // auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector> >();
    // subproblem->pc_type("asm");
    // TrustRegion<Matrix, Vector> solver(subproblem);
    // subproblem->atol(1e-14);
    // solver.verbose(true);
    // solver.delta0(1e-0);
    // solution.set(1.0);
    // solver.atol(1e-10);
    // solver.solve(*level_functions_.back(), solution);

    // auto subproblem = std::make_shared<SteihaugToint<Matrix, Vector> >();
    // subproblem->pc_type("asm");
    // Newton<Matrix, Vector> solver(subproblem);
    // solver.verbose(true);
    // solution.set(1.0);
    // solver.solve(*level_functions_.back(), solution);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    write_to_file(*spaces_[n_levels_ - 1], solution);

  }

 private:
  bool init_;
  SizeType n_levels_;
  SizeType n_coarse_sub_comm_;

  std::vector<std::shared_ptr<FunctionSpace>> spaces_;
  std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> transfers_;

  std::vector<std::shared_ptr<ExtendedFunction<Matrix, Vector>>>
      level_functions_;
  std::vector<std::shared_ptr<BCType>> BC_conditions_;


  std::string log_output_path_;
  std::string output_path_;

  std::shared_ptr<
      RMTR_inf<Matrix, Vector, TRBoundsGratton<Matrix, Vector>, SECOND_ORDER> >
      rmtr_;

  bool save_output_;
};

}  // namespace utopia

#endif  // UTOPIA_DM_RMTR_SETUP_STEADY_STATE_HPP

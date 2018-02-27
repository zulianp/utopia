#ifndef UTOPIA_PETSC_MULTIGRID_HPP
#define UTOPIA_PETSC_MULTIGRID_HPP 

#include "utopia_Multigrid.hpp"
#include "utopia_IterativeSolver.hpp"

namespace utopia {

	template<class Matrix, class Vector>
	class Multigrid<Matrix, Vector, PETSC_EXPERIMENTAL> : public IterativeSolver<Matrix, Vector>, public MultiLevelBase<Matrix, Vector> {
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

		typedef utopia::LinearSolver<Matrix, Vector>        Solver;
		typedef utopia::IterativeSolver<Matrix, Vector>     IterativeSolver;
		typedef utopia::Smoother<Matrix, Vector>            Smoother;
		typedef utopia::Level<Matrix, Vector>               Level;
		typedef utopia::Transfer<Matrix, Vector>            Transfer;
	public:
		void update(const std::shared_ptr<const Matrix> &op) override
		{
			if(!ksp_) {
				init_ksp(op);
			}

			this->galerkin_assembly(op);
			KSPSetOperators(*ksp_, raw_type(*op), raw_type(*op));

			PC pc;
			KSPGetPC(*ksp_, &pc);
			for (std::size_t i = 0; i < this->num_levels()-1; i++)
			{
				KSP smoother;
				PCMGGetSmoother(pc, i, &smoother);
				KSPSetOperators(smoother, raw_type(this->level(i).A()), raw_type(this->level(i).A()));
			}
		}

		virtual bool apply(const Vector &rhs, Vector &x) override
		{
			if(this->verbose())
				this->init_solver("utopia/petsc Multigrid",  {" it.", "|| Au - b||"});

			KSPSolve(*ksp_, raw_type(rhs), raw_type(x));
			
			KSPConvergedReason  reason;
			PetscInt            its; 
			KSPGetConvergedReason(*ksp_, &reason);
			KSPGetIterationNumber(*ksp_, &its);
			
			if(this->verbose())
				this->exit_solver(its, reason); 
			//FIXME
			return true;
		}

		Multigrid(const std::shared_ptr<Smoother> &smoother    = nullptr,
		          const std::shared_ptr<Solver> &linear_solver = nullptr,
		          const Parameters params = Parameters())
		: smoother_(smoother), linear_solver_(linear_solver) 
		{
		    set_parameters(params); 
		}

		void set_parameters(const Parameters params) override
		{
		    IterativeSolver::set_parameters(params); 
		    MultiLevelBase<Matrix, Vector>::set_parameters(params); 
		}

	private:
		std::shared_ptr<Smoother> smoother_;
		std::shared_ptr<Solver>   linear_solver_;

		std::shared_ptr<KSP> ksp_;
		std::vector<Level> levels_;

		template<class Solver>
		bool set_solver(Solver &solver, KSP &ksp)
		{
			// auto * casted = dynamic_cast< KSPSolver<Matrix, Vector> *>(&solver);

			// if(casted) {
			// 	KSPSetType(ksp, casted->ksp_type().c_str());
			// 	KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);

			// 	PC pc; 
			// 	KSPGetPC(ksp, &pc);

			// 	if(!casted->get_preconditioner()) {
			// 	    PCSetType(pc, casted->pc_type().c_str());
			// 	} else {
			// 		PCSetType(pc, PCSOR);
			// 	}

			// 	return true;

			// } else {
			// 	// auto * factor = dynamic_cast< Factorization<Matrix, Vector> *>(&solver);
				
			// 	// if(factor) {
			// 	// 	factor->strategy().set_ksp_options(ksp);
			// 	// 	KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
			// 	// 	return true;
			// 	// } 

			// 	// else {
			// 	// 	return false;
			// 	// }
			// }

			return false;
		}

		void init_ksp(const std::shared_ptr<const Matrix> &op)
		{
			auto comm = op->implementation().communicator();
			ksp_ = std::shared_ptr<KSP>(new KSP, [](KSP *ksp) { KSPDestroy(ksp); delete ksp; ksp = nullptr; });

			KSPCreate(comm, ksp_.get());
			KSPSetFromOptions(*ksp_);

			PC pc;
			KSPGetPC(*ksp_, &pc);
			PCSetType(pc, PCMG);

			PCMGSetLevels(pc, this->num_levels(), nullptr);
			// PCMGSetGalerkin(pc, PETSC_TRUE);
			PCMGSetGalerkin(pc, PETSC_FALSE);
			KSPSetInitialGuessNonzero(*ksp_, PETSC_TRUE);

			for (std::size_t i = 0; i < this->num_levels()-1; i++)
			{
				KSP smoother;
				PC sm;
				PCMGGetSmoother(pc, i + 1, &smoother);

				bool user_solver = false;

				if(i > 0 && smoother_) {
					user_solver = set_solver(*smoother_, smoother);
				} else if(linear_solver_) {
					user_solver = set_solver(*linear_solver_, smoother);
				} 

				if(!user_solver) {
					std::cerr << "[Warning] not using user solvers" << std::endl;
					
					KSPSetType(smoother, KSPRICHARDSON);
					KSPGetPC(smoother, &sm);
					PCSetType(sm, PCSOR);

					// KSPSetInitialGuessNonzero(smoother, PETSC_TRUE);
				}

				Mat I = raw_type(this->transfer(i).I());
				PCMGSetInterpolation(pc, i + 1, I);
			}

			PCMGSetNumberSmoothUp(pc,   this->post_smoothing_steps());
			PCMGSetNumberSmoothDown(pc, this->pre_smoothing_steps());

			KSPSetTolerances(*ksp_, this->rtol(), this->atol(), PETSC_DEFAULT, this->max_it());

			if(this->verbose()) {
				KSPMonitorSet(
				*ksp_,
				[](KSP, PetscInt iter, PetscReal res, void*) -> PetscErrorCode{
					PrintInfo::print_iter_status({static_cast<Scalar>(iter), res}); 
					return 0;
				},
				nullptr,
				nullptr);
			}
		}
	};

}

#endif //UTOPIA_PETSC_MULTIGRID_HPP

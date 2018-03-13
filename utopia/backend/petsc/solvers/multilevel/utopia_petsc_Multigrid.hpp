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
		: smoother_(smoother), linear_solver_(linear_solver), default_ksp_type_(KSPRICHARDSON), default_pc_type_(PCSOR)
		{
		    set_parameters(params); 
		}

		void set_parameters(const Parameters params) override
		{
		    IterativeSolver::set_parameters(params); 
		    MultiLevelBase<Matrix, Vector>::set_parameters(params); 
		}

		inline void set_default_ksp_type(const KSPType &ksp_type)
		{
			default_ksp_type_ = ksp_type;
		}

		inline void set_default_pc_type(const PCType &pc_type)
		{
			default_pc_type_ = pc_type;
		}

	private:
		std::shared_ptr<Smoother> smoother_;
		std::shared_ptr<Solver>   linear_solver_;

		std::shared_ptr<KSP> ksp_;
		std::vector<Level> levels_;

		KSPType default_ksp_type_;
		PCType default_pc_type_;

		template<class Solver>
		bool set_solver(Solver &solver, KSP &ksp)
		{
			auto * casted = dynamic_cast< KSPSolver<Matrix, Vector> *>(&solver);

			if(casted) {
				// if(casted->ksp_type() == KSPCG || casted->ksp_type() == KSPFCG) {
				// 	KSPSetType(ksp, KSPFCG);	
				// } else 
				if(casted->ksp_type() == KSPGMRES || casted->ksp_type() == KSPFGMRES) {
					KSPSetType(ksp, KSPFGMRES);	
				} else {
					if(casted->ksp_type() != default_ksp_type_) {
						std::cerr << "[Warning] " << casted->ksp_type() << " not supported by petsc PCMG using fallback option KSPRICHARDSON" << std::endl;
					}

					KSPSetType(ksp, default_ksp_type_);
				}
				
				PC pc; 
				KSPGetPC(ksp, &pc);

				if(!casted->get_preconditioner()) {
				    PCSetType(pc, casted->pc_type().c_str());
				} else {
				// 	//FIXME use PCShell
					PCSetType(pc, default_pc_type_);
				}

				{
					PCType pc_type;
					PCGetType(pc, &pc_type);

					KSPType ksp_type;
					KSPGetType(ksp, &ksp_type);

					std::cout << ksp_type << ", " << pc_type << std::endl;
				}

				// return false;
				return true;

			} else {
				auto * factor = dynamic_cast< Factorization<Matrix, Vector, PETSC> *>(&solver);
				
				if(factor) {
					factor->strategy().set_ksp_options(ksp);
					KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
					return true;
				} 

				// else {
				// 	return false;
				// }
			}

			return false;
		}

		void init_ksp(const std::shared_ptr<const Matrix> &op)
		{
			auto comm = op->implementation().communicator();
			ksp_ = std::shared_ptr<KSP>(new KSP, [](KSP *&ksp) { KSPDestroy(ksp); delete ksp; ksp = nullptr; });

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

				if(smoother_) {
					user_solver = set_solver(*smoother_, smoother);
					if(!user_solver) {
						std::cerr << "[Warning] not using user smoother" << std::endl;
					}
				} 

				if(!user_solver) {
					// KSPFGMRES, KSPGCG, or KSPRICHARDSON 
					KSPSetType(smoother, default_ksp_type_);
					KSPGetPC(smoother, &sm);
					PCSetType(sm, default_pc_type_);

					// KSPSetInitialGuessNonzero(smoother, PETSC_TRUE);
				}

				Mat I = raw_type(this->transfer(i).I());
				PCMGSetInterpolation(pc, i + 1, I);

				if(linear_solver_) {
					KSP coarse_solver;
					PCMGGetCoarseSolve(pc, &coarse_solver);
					bool user_solver = set_solver(*linear_solver_, coarse_solver);

					if(!user_solver) {
						std::cerr << "[Warning] not using user linear_solver" << std::endl;
					}
				}
			}

			PCMGSetNumberSmoothUp(pc,   this->post_smoothing_steps());
			PCMGSetNumberSmoothDown(pc, this->pre_smoothing_steps());

			KSPSetTolerances(*ksp_, this->rtol(), this->atol(), PETSC_DEFAULT, this->max_it());

			if(this->verbose()) {
				KSPMonitorSet(
				*ksp_,
				[](KSP, PetscInt iter, PetscReal res, void*) -> PetscErrorCode {
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

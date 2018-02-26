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
			this->init_solver("utopia/petsc Multigrid",  {" it.", "|| Au - b||"});

			KSPSolve(*ksp_, raw_type(rhs), raw_type(x));
			
			KSPConvergedReason  reason;
			PetscInt            its; 
			KSPGetConvergedReason(*ksp_, &reason);
			KSPGetIterationNumber(*ksp_, &its);
			
			this->exit_solver(its, reason); 
			//FIXME
			return true;
		}

		Multigrid(const std::shared_ptr<Smoother> &smoother    = nullptr,
		          const std::shared_ptr<Solver> &direct_solver = nullptr,
		          const Parameters params = Parameters())
		//: _smoother(smoother), _direct_solver(direct_solver) 
		{
			std::cerr << "[Warning] this Multigrid is not using user provided smoother and direct solver" << std::endl;
		    set_parameters(params); 
		}

		void set_parameters(const Parameters params) override
		{
		    IterativeSolver::set_parameters(params); 
		    MultiLevelBase<Matrix, Vector>::set_parameters(params); 
		}

	private:
		std::shared_ptr<KSP> ksp_;
		std::vector<Level> levels_;

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

			for (std::size_t i = 0; i < this->num_levels()-1; i++)
			{
				KSP smoother;
				PC sm;
				PCMGGetSmoother(pc, i + 1, &smoother);
				KSPSetType(smoother, KSPRICHARDSON);
				KSPGetPC(smoother, &sm);
				PCSetType(sm, PCSOR);

				Mat I = raw_type(this->transfer(i).I());
				PCMGSetInterpolation(pc, i + 1, I);
			}

			PCMGSetNumberSmoothUp(pc, 3);
			PCMGSetNumberSmoothDown(pc, 3);


			KSPSetInitialGuessNonzero(*ksp_, PETSC_TRUE);
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

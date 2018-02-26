#ifndef UTOPIA_PETSC_MULTIGRID_HPP
#define UTOPIA_PETSC_MULTIGRID_HPP 

#include "utopia_Multigrid.hpp"
#include "utopia_IterativeSolver.hpp"

namespace utopia {

	template<class Matrix, class Vector>
	class Multigrid<Matrix, Vector, PETSC_EXPERIMENTAL> : public IterativeSolver<Matrix, Vector> {
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

			KSPSetOperators(*ksp_, raw_type(*op), raw_type(*op));
		}

		virtual bool apply(const Vector &rhs, Vector &x) override
		{
			KSPSolve(*ksp_, raw_type(rhs), raw_type(x));
			
			KSPConvergedReason  reason;
			PetscInt            its; 
			KSPGetConvergedReason(*ksp_, &reason);
			KSPGetIterationNumber(*ksp_, &its);
			this->exit_solver(its, reason); 

			//FIXME
			return true;
		}

		virtual bool init_interpolators_from_coarse_to_fine(const std::vector<std::shared_ptr<Matrix>> &interpolation)
		{
			transfers_.clear();
			transfers_.reserve(interpolation.size());
			
			for(auto I_ptr : interpolation) {
				transfers_.push_back(Transfer(I_ptr));
			}

			return !transfers_.empty();
		}

		

		inline std::size_t num_levels() const
		{
			return transfers_.size() + 1;
		}

	private:
		std::shared_ptr<KSP> ksp_;
		std::vector<Transfer> transfers_;

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
			PCMGSetGalerkin(pc, PETSC_TRUE);

			for (std::size_t i = 0; i < transfers_.size(); i++)
			{
				KSP smoother;
				PC sm;
				PCMGGetSmoother(pc, i + 1, &smoother);
				KSPSetType(smoother, KSPRICHARDSON);
				KSPGetPC(smoother, &sm);
				PCSetType(sm, PCSOR);

				Mat I = raw_type(transfers_[i].I());
				PCMGSetInterpolation(pc, i + 1, I);
			}

			PCMGSetNumberSmoothUp(pc, 3);
			PCMGSetNumberSmoothDown(pc, 3);
		}
	};
}

#endif //UTOPIA_PETSC_MULTIGRID_HPP

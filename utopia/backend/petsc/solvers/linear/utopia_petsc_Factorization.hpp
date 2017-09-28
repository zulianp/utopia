#ifndef UTOPIA_FACTORIZATION_HPP
#define UTOPIA_FACTORIZATION_HPP  

#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_DirectSolver.hpp"

namespace utopia {

	//FIXME
	enum DirectSolverLib {
		PETSC_TAG = 0,

#ifdef PETSC_HAVE_MUMPS
		MUMPS_TAG = 1,
#endif //WITH_MUMPS   

#ifdef PETSC_HAVE_SUPERLU
		SUPERLU_TAG = 2,
#endif //WITH_SUPERLU

#ifdef PETSC_HAVE_SUPERLU_DIST
		SUPERLU_DIST_TAG = 3,
#endif //WITH_SUPERLU_DIST    	

	};
	
	static const int N_SOLVER_LIBS = 4;

//////////////////////////////////////////////////////////

	enum DirectSolverType {
		LU_DECOMPOSITION_TAG = 0,
		CHOLESKY_DECOMPOSITION_TAG = 1
	};
	
	static const int N_SOLVER_TYPES = 2;

//////////////////////////////////////////////////////////

	static const char * SOLVER_LIBS [N_SOLVER_LIBS] = { "petsc", "mumps", "superlu", "superlu_dist"};
	static const char * SOLVER_TYPES[N_SOLVER_TYPES]   = {"lu", "cholesky" };  


//////////////////////////////////////////////////////////
	template<typename Matrix, typename Vector>
	class Factorization<Matrix, Vector, PETSC> : public DirectSolver<Matrix, Vector>
	{
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;

	public:
		Factorization(const std::string sp = "mumps", const std::string pct = "lu", const Parameters params = Parameters()):
		strategy_(params,   {" "}, 
			{"lu", "jacobi", "sor", "shell",  "bjacobi",  "ilu",  "icc", "cholesky", "pbjacobi"}, 
			{"mumps", "superlu", "superlu_dist", "petsc", "cusparse"} )
		{ 
			if(pct=="lu" && sp == "mumps")
				set_type(MUMPS_TAG, LU_DECOMPOSITION_TAG);
			else if(pct=="cholesky")
				set_type(PETSC_TAG, CHOLESKY_DECOMPOSITION_TAG);
			else
				set_type(PETSC_TAG, LU_DECOMPOSITION_TAG);
		}

		void set_type(DirectSolverLib lib, DirectSolverType type)
		{
			strategy_.solver_package(SOLVER_LIBS[lib]); 
			strategy_.pc_type(SOLVER_TYPES[type]); 
		}

		inline bool apply(const Vector &b, Vector &x) override
		{
			return strategy_.apply(b, x);
		}  

		inline void update(const std::shared_ptr<const Matrix> &op) override
		{
			strategy_.update(op);
		}     

		inline void set_parameters(const Parameters params) override
		{
			LinearSolver<Matrix, Vector>::set_parameters(params);
			strategy_.set_parameters(params);
		}                                           

	private:

		class Strategy : public KSPSolver<Matrix, Vector> {
		public:
			using KSPSolver<Matrix, Vector>::KSPSolver;
	            /*@todo use : PCFactorSetReuseFill(PC pc,PetscBool flag)  // to keep factorization from prev. levels 
	            */

			void set_ksp_options(KSP & ksp) override
			{
				this->reset_preconditioner(); 
				PetscErrorCode ierr;
				ierr = KSPSetType(ksp, KSPPREONLY);
				PC pc; 
				ierr = KSPGetPC(ksp,&pc);

				ierr = PCSetType(pc, this->pc_type().c_str());
				ierr = PCFactorSetMatSolverPackage(pc, this->solver_package().c_str());
				ierr = KSPSetTolerances(ksp, IterativeSolver::rtol(), IterativeSolver::atol(), PETSC_DEFAULT,  IterativeSolver::max_it());
			}
		};

		Strategy strategy_;
	};
}

#endif //UTOPIA_FACTORIZATION_HPP

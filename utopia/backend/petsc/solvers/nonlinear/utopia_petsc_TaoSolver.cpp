#include "utopia_petsc_TaoSolver.hpp"
#include "petsctao.h"
#include "utopia_petsc_Vector.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Types.hpp"

#define U_CHECKERR(ierr) { if(ierr != 0) return false; }

namespace utopia {

	template class TaoSolver<DSMatrixd, DVectord>;
	template class TaoSolver<DMatrixd, DVectord>;

	template<class Matrix, class Vector>
	static PetscErrorCode UtopiaTaoEvaluateObjective(Tao tao, Vec x, PetscReal *ret, void *ctx)
	{
		using Scalar = typename Function<Matrix, Vector>::Scalar;

		Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);
		assert(fun);

		Vector utopia_x;
		convert(x, utopia_x);

		Scalar utopia_value = 0.;
		if(!fun->value(utopia_x, utopia_value)) {
			return 1;
		}

		*ret = utopia_value;
		return 0;
	}

	template<class Matrix, class Vector>
	static PetscErrorCode UtopiaTaoEvaluateGradient(Tao tao, Vec x, Vec g, void *ctx)
	{ 
		Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);
		assert(fun);

		Vector utopia_x;
		convert(x, utopia_x);

		Vector utopia_g = local_zeros(local_size(utopia_x));
		if(!fun->gradient(utopia_x, utopia_g)) {
			return 1;
		}

		convert(utopia_g, g);
		return 0;
	}

	template<class Matrix, class Vector>
	static PetscErrorCode UtopiaTaoFormHessian(Tao tao, Vec x, Mat H, Mat Hpre, void *ctx)
	{
		// PetscErrorCode ierr = 0;
		Function<Matrix, Vector> * fun = static_cast<Function<Matrix, Vector> *>(ctx);
		assert(fun);

		Vector utopia_x;

		Matrix utopia_H;
		Matrix utopia_Hpre;
		
		convert(x, utopia_x);
		utopia_H.implementation().wrap(H);
		utopia_Hpre.implementation().wrap(Hpre);

		if(!fun->hessian(utopia_x, utopia_H, utopia_Hpre)) {
			if(!fun->hessian(utopia_x, utopia_H)) {
				utopia_error("[Error] Failed to assemble Hessian.");
				assert(false);
				return 1;
			} 

		} else {
			if(raw_type(utopia_Hpre) != Hpre) {
				//FIXME maybe add optimization options
				MatCopy(raw_type(utopia_Hpre), Hpre, DIFFERENT_NONZERO_PATTERN);
			}
		}

		if(raw_type(utopia_H) != H) {
			//FIXME maybe add optimization options
			MatCopy(raw_type(utopia_H), H, DIFFERENT_NONZERO_PATTERN);
		}

		return 0;
	}

	template<class Matrix, class Vector>
	bool UtopiaTaoSetUp(Tao tao, Function<Matrix, Vector> &fun)
	{
		fun.data()->init();
		if(!fun.initialize_hessian(*fun.data()->H, *fun.data()->H_pre)) {
			utopia_error("TaoSolver requires Function::initialize_hessian to be implemented.");
			assert(false);
			return false;
		}

		PetscErrorCode ierr = 0;
		if(fun.has_preconditioner()) 
		{
			ierr = TaoSetHessianRoutine(
				tao,
				raw_type(*fun.data()->H),
				raw_type(*fun.data()->H_pre),
				UtopiaTaoFormHessian<Matrix, Vector>,
			    static_cast<void *>(&fun)); U_CHECKERR(ierr);

		} else {
			ierr = TaoSetHessianRoutine(
				tao,
				raw_type(*fun.data()->H),
				raw_type(*fun.data()->H),
				UtopiaTaoFormHessian<Matrix, Vector>,
			    static_cast<void *>(&fun)); U_CHECKERR(ierr);
		}

		ierr = TaoSetObjectiveRoutine(tao,
			UtopiaTaoEvaluateObjective<Matrix, Vector>,
			static_cast<void *>(&fun)); U_CHECKERR(ierr);

		ierr = TaoSetGradientRoutine(
			tao,
			UtopiaTaoEvaluateGradient<Matrix, Vector>,
			static_cast<void *>(&fun)); U_CHECKERR(ierr);

		return false;
	} 

	void TaoSolverWrapper::set_function(Function<DSMatrixd, DVectord> &fun)
	{
		auto tao = (Tao *) &data_;
		UtopiaTaoSetUp(*tao, fun);
	}

	void TaoSolverWrapper::set_function(Function<DMatrixd, DVectord> &fun)
	{
		auto tao = (Tao *) &data_;
		UtopiaTaoSetUp(*tao, fun);
	}

	bool TaoSolverWrapper::get_ksp(KSP *ksp)
	{
		if(!data_) return false;

		auto tao = (Tao *) &data_;
		PetscErrorCode ierr = TaoGetKSP(*tao, ksp); U_CHECKERR(ierr);
		return true;
	}

	bool TaoSolverWrapper::init(
		MPI_Comm comm,
		const std::string &type,
		const PetscReal gatol,
		const PetscReal grtol,
		const PetscReal gttol,
		const PetscInt maxits)
	{
		auto tao = (Tao *) &data_;
		PetscErrorCode ierr = 0;
		ierr = TaoCreate(comm, tao);   U_CHECKERR(ierr);
		
		if(type.empty()) {
			ierr = TaoSetType(*tao, TAOTRON); U_CHECKERR(ierr);
		} else {
			ierr = TaoSetType(*tao, type.c_str()); U_CHECKERR(ierr);
		}

		ierr = TaoSetTolerances(*tao, gatol, grtol, gttol); U_CHECKERR(ierr);
		ierr = TaoSetMaximumIterations(*tao, maxits); U_CHECKERR(ierr);

		KSP ksp;
		PC pc;

		ierr = TaoGetKSP(*tao, &ksp); U_CHECKERR(ierr);

		if(ksp) {
			ierr = KSPSetType(ksp, ksp_type_.c_str()); U_CHECKERR(ierr);

			ierr = KSPGetPC(ksp, &pc); U_CHECKERR(ierr);
			ierr = PCSetType(pc, pc_type_.c_str()); U_CHECKERR(ierr);
	
#if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
			ierr = PCFactorSetMatSolverPackage(pc, solver_package_.c_str()); U_CHECKERR(ierr);
#else
			m_utopia_error("PCFactorSetMatSolverPackage not available in petsc 3.9.0 find equivalent");
#endif 
			ierr = KSPSetInitialGuessNonzero(ksp, PETSC_FALSE); U_CHECKERR(ierr);
		} else {
			utopia_error("Tao does not have a ksp");
			assert(false);
			return false;
		}

		ierr = TaoSetFromOptions(*tao);  U_CHECKERR(ierr);
		return true;
	}

	void TaoSolverWrapper::destroy()
	{
		if(data_) {
			Tao * tao = (Tao *) &data_;
			TaoDestroy(tao);
		}
	}

	void TaoSolverWrapper::set_monitor(MPI_Comm comm)
	{
		const char  monfilename[7] ="stdout";
		PetscViewer    monviewer;
		auto tao = (Tao *) &data_;
		PetscViewerASCIIOpen(comm, monfilename, &monviewer);
		TaoSetMonitor(*tao, TaoDefaultSMonitor, monviewer, (PetscErrorCode (*)(void**))PetscViewerDestroy);
	}


	TaoSolverWrapper::TaoSolverWrapper()
	: data_(nullptr), ksp_type_(KSPPREONLY), pc_type_(PCLU), solver_package_("mumps")
	{}

	TaoSolverWrapper::~TaoSolverWrapper()
	{
		destroy();
	}

	void TaoSolverWrapper::set_ksp_types(const std::string &ksp, const std::string &pc, const std::string &solver_package)
	{
		ksp_type_ = ksp;
		pc_type_ = pc;
		solver_package_ = solver_package;
	}

	void TaoSolverWrapper::set_pc_type(const std::string &pc)
	{
		pc_type_ = pc;
	}


	bool TaoSolverWrapper::set_bounds(const PetscVector &lb, const PetscVector &ub)
	{
		auto tao = static_cast<Tao>(data_);

		PetscErrorCode ierr = 0; 
		// ierr = TaoSetInequalityBounds(tao, lb.implementation(), ub.implementation()); U_CHECKERR(ierr);
		ierr = TaoSetVariableBounds(tao, lb.implementation(), ub.implementation()); U_CHECKERR(ierr);
		return true;
	}

	bool TaoSolverWrapper::solve(PetscVector &x)
	{
		auto tao = static_cast<Tao>(data_);

		PetscErrorCode ierr = 0; 
		TaoSetInitialVector(tao, x.implementation());
		ierr = TaoSolve(tao); U_CHECKERR(ierr);

		PetscInt iterate;
		PetscReal f;
		PetscReal gnorm;
		PetscReal cnorm;
		PetscReal xdiff;
		TaoConvergedReason reason;
		TaoGetSolutionStatus(tao, &iterate, &f, &gnorm, &cnorm, &xdiff, &reason);

		// if(this->verbose()) {
			// std::cout << "iterate: " << iterate << std::endl;
			// std::cout << "f: " << f << std::endl;
			// std::cout << "gnorm: " << gnorm << std::endl;
			// std::cout << "cnorm: " << cnorm << std::endl;
			// std::cout << "xdiff: " << xdiff << std::endl;
			// std::cout << "reason: " << reason << std::endl;
		// }

		if(reason < 0) {
			utopia_error("> Failed to converge");
		}

		return reason >= 0;
	}
}

#undef U_CHECKERR

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
		std::cout << "UtopiaTaoFormHessian" << std::endl;
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
				std::cerr << "HERERE" << std::endl;
				return 1;
			} 
			// else {
			// 	PetscBool assembled = PETSC_FALSE;
			// 	PetscErrorCode ierr = MatAssembled(Hpre, &assembled); CHKERRQ(ierr);

			// 	if(!assembled) {
			// 		std::cerr << "[Error] Handle outside!" << std::endl;
			// 	}

			// 	MatCopy(H, Hpre, SAME_NONZERO_PATTERN); 
			// }
		}
#ifndef NDEBUG
		else {
			assert(raw_type(utopia_Hpre) == Hpre);
		}
#endif


		assert(raw_type(utopia_H) == H);
		return 0;
	}

	template<class Matrix, class Vector>
	bool UtopiaTaoSetUp(Tao tao, Function<Matrix, Vector> &fun)
	{
		PetscErrorCode ierr = 0;
		if(fun.has_preconditioner()) 
		{
			fun.data()->init();
			ierr = TaoSetHessianRoutine(
				tao,
				raw_type(*fun.data()->H),
				raw_type(*fun.data()->H_pre),
				UtopiaTaoFormHessian<Matrix, Vector>,
			    static_cast<void *>(&fun)); U_CHECKERR(ierr);

		} else {
			fun.data()->init();
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

	// template<>
	// void TaoSolverWrapper::set_function<DMatrixd, DVectord>(Function<DSMatrixd, DVectord> &fun)
	// {
	// 	auto tao = (Tao *) &data_;
	// 	UtopiaTaoSetUp(tao, fun);
	// }

	// UtopiaTaoSetHessianRoutine(
	// 	Tao tao,
	// 	Mat H,
	// 	Mat Hpre,
	// 	PetscErrorCode (*FormHessian)(Tao,Vec,Mat,Mat,
	// void*), void *user);

	bool TaoSolverWrapper::init(MPI_Comm comm)
	{
		auto tao = (Tao *) &data_;
		PetscErrorCode ierr = 0;
		ierr = TaoCreate(comm, tao);    U_CHECKERR(ierr);
		ierr = TaoSetFromOptions(*tao); U_CHECKERR(ierr);
		return true;
	}

	void TaoSolverWrapper::destroy()
	{
		if(data_) {
			Tao * tao = (Tao *) &data_;
			TaoDestroy(tao);
		}
	}

	TaoSolverWrapper::TaoSolverWrapper()
	: data_(nullptr)
	{}

	TaoSolverWrapper::~TaoSolverWrapper()
	{
		destroy();
	}

	bool TaoSolverWrapper::set_bounds(const PetscVector &lb, const PetscVector &ub)
	{
		auto tao = static_cast<Tao>(data_);

		PetscErrorCode ierr = 0; 
		ierr = TaoSetInequalityBounds(tao, lb.implementation(), ub.implementation()); U_CHECKERR(ierr);

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
		// 	std::cout << "iterate: " << iterate << std::endl;
		// 	std::cout << "f: " << f << std::endl;
		// 	std::cout << "gnorm: " << gnorm << std::endl;
		// 	std::cout << "cnorm: " << cnorm << std::endl;
		// 	std::cout << "xdiff: " << xdiff << std::endl;
		// 	std::cout << "reason: " << reason << std::endl;
		// }
		
		return reason >= 0;
	}
}

#undef U_CHECKERR
#ifndef UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP
#define UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP

#include "utopia_SemismoothNewton.hpp"
#include "utopia_PETScKSPSolver.hpp"

// PetscErrorCode  KSPRegister(const char sname[],PetscErrorCode (*function)(KSP))
// KSPRegister("my_solver",MySolverCreate);

namespace utopia {

	//FIXME and then add PETSC to the backend flag
	static const int PETSC_EXPERIMENTAL = -1000;

	template<class Matrix, class Vector>
	class SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL> : public IterativeSolver<Matrix, Vector> {
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::LinearSolver<Matrix, Vector> Solver;
		typedef utopia::BoxConstraints<Vector>      BoxConstraints;
		
	public:
		
		SemismoothNewton(const std::shared_ptr <Solver> &linear_solver   = std::shared_ptr<Solver>(),
						 const Parameters params                         = Parameters() ) :
		linear_solver_(linear_solver)
		{
			set_parameters(params);
		}
				
		bool solve(const Matrix &A, const Vector &b, Vector &x)  override
		{
			PetscErrorCode ierr = 0.;
			Vector f = local_zeros(local_size(b));
			Matrix J = A;

			SemismoothNewtonCtx ctx;
			ctx.H = &A;
			ctx.g = &b;

			Vec lobo, upbo;

			DVectord dummy_lobo, dummy_upbo;
			if(constraints_.has_upper_bound()) {
				upbo = raw_type(*constraints_.upper_bound());
			} else {
				dummy_upbo = local_values(local_size(b).get(0), PETSC_INFINITY);
				upbo = raw_type(dummy_upbo);
			}

			if(constraints_.has_lower_bound()) {
				lobo = raw_type(*constraints_.lower_bound());
			}  else {
				dummy_lobo = local_values(local_size(b).get(0), PETSC_NINFINITY);
				lobo = raw_type(dummy_lobo);
			}

			linear_solver_->update(make_ref(A));

			SNES snes;
			SNESCreate(PETSC_COMM_WORLD, &snes);
			SNESSetType(snes, SNESVINEWTONSSLS);
			SNESSetFunction(snes, raw_type(f), SemismoothNewton::Gradient, &ctx);
			SNESSetJacobian(snes, raw_type(J), raw_type(J), SemismoothNewton::Hessian, &ctx);
			SNESVISetVariableBounds(snes, lobo, upbo);

			KSP ksp;
			PC pc; 
			SNESGetKSP(snes, &ksp);
			KSPSetType(ksp, KSPPREONLY);		
			KSPGetPC(ksp, &pc);

			PCSetType(pc, "lu");
			PCFactorSetMatSolverPackage(pc, "mumps");


			//FIXME use utopia linear solvers
			// auto shell_ptr = linear_solver_.get();
			// assert(shell_ptr);
			// ierr = PCSetType(pc, PCSHELL);  
			// ierr = PCShellSetApply(pc, UtopiaPCApplyShell);
			// ierr = PCShellSetContext(pc, shell_ptr);
			// ierr = PCShellSetName(pc, "Utopia Linear Solver");

			//
			ierr = KSPSetTolerances(ksp, this->rtol(), this->atol(), PETSC_DEFAULT, this->max_it());

			SNESSetFromOptions(snes);
			SNESSolve(snes, nullptr, raw_type(x));

			SNESDestroy(&snes);
			return true;
		}
		
		virtual void set_parameters(const Parameters params) override
		{
			IterativeSolver<Matrix, Vector>::set_parameters(params);
		}
		
		virtual bool set_box_constraints(const BoxConstraints & box)
		{
			constraints_ = box;
			return true;
		}
		
		virtual bool get_box_constraints(BoxConstraints & box)
		{
			box = constraints_;
			return true;
		}
		
	private:	
		std::shared_ptr <Solver>        linear_solver_;
		BoxConstraints                  constraints_;

		typedef struct {
			const Matrix * H;
			const Vector  * g;
		} SemismoothNewtonCtx; 

		static PetscErrorCode Gradient(SNES snes, Vec x, Vec f, void*ctx)
		{
			SemismoothNewtonCtx * ssn_ctx = (SemismoothNewtonCtx *) ctx;
			Vector x_utopia, f_utopia;
			convert(x, x_utopia);
			f_utopia  = (*ssn_ctx->H) * x_utopia - (*ssn_ctx->g);
			convert(f_utopia, f);
			return 0;
		}

		static PetscErrorCode Hessian(SNES snes, Vec x, Mat Hmat, Mat Pmat, void*ctx)
		{
			SemismoothNewtonCtx * ssn_ctx = (SemismoothNewtonCtx *) ctx;
			assert(Hmat == Pmat);
			convert((*ssn_ctx->H), Hmat);
			return 0;
		}
	};

}

#endif //UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP

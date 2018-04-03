#ifndef UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP
#define UTOPIA_PETSC_SEMI_SMOOTH_NEWTON_HPP

#include "utopia_SemismoothNewton.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include <petscsnes.h>
#include "utopia_petsc.hpp"

// PetscErrorCode  KSPRegister(const char sname[],PetscErrorCode (*function)(KSP))
// KSPRegister("my_solver",MySolverCreate);

namespace utopia {

	template<class Matrix, class Vector>
	class SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL> : public IterativeSolver<Matrix, Vector> {
		typedef UTOPIA_SCALAR(Vector)    Scalar;
		typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
		typedef utopia::LinearSolver<Matrix, Vector> Solver;
		typedef utopia::BoxConstraints<Vector>      BoxConstraints;
		
	public:
		
		SemismoothNewton(const std::shared_ptr <Solver> &linear_solver   = std::shared_ptr<Solver>(),
						 const Parameters params                         = Parameters() ) :
		linear_solver_(linear_solver), line_search_type_(SNESLINESEARCHBASIC)
		{
			set_parameters(params);
		}
				
		bool solve(const Matrix &A, const Vector &b, Vector &x)  override
		{
			Vector f = local_zeros(local_size(b));
			Matrix J = A;

			if(empty(x)) {
				x = local_zeros(local_size(b));
			}

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

			SNES snes;
			SNESCreate(A.implementation().communicator(), &snes);
			SNESSetFromOptions(snes);
			SNESSetType(snes, SNESVINEWTONSSLS);
			SNESSetFunction(snes, raw_type(f), SemismoothNewton::Gradient, &ctx);
			SNESSetJacobian(snes, raw_type(J), raw_type(J), SemismoothNewton::Hessian, &ctx);
			SNESVISetVariableBounds(snes, lobo, upbo);
#if !UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)  
			// asdf
			SNESSetForceIteration(snes, PETSC_TRUE);
#endif

			KSP ksp;
			PC pc; 

			SNESGetKSP(snes, &ksp);
			KSPGetPC(ksp, &pc);

			auto ksp_solver_ptr = std::dynamic_pointer_cast<KSPSolver<DSMatrixd, DVectord> >(linear_solver_);
			
			bool has_linear_solver = false;
			
			if(ksp_solver_ptr) {
				ksp_solver_ptr->set_ksp_options(ksp);

				if(ksp_solver_ptr->verbose()) {
				    KSPMonitorSet(ksp,
								[](KSP, PetscInt iter, PetscReal res, void*) -> PetscErrorCode {
									PrintInfo::print_iter_status({static_cast<Scalar>(iter), res}); 
									return 0;
								},
								nullptr,
								nullptr);
				}

				has_linear_solver = true;

			} else {
				auto factor_solver_ptr = std::dynamic_pointer_cast< Factorization<Matrix, Vector, PETSC> >(linear_solver_);
				if(factor_solver_ptr) {
					factor_solver_ptr->strategy().set_ksp_options(ksp);
					has_linear_solver = true;

					KSPSetTolerances(ksp, 0, 0, PETSC_DEFAULT, 1);
				}
			}

			if(!has_linear_solver) {
				std::cerr << "Non-petsc linear solvers not supported yet: falling-back to mumps/lu" << std::endl;
				KSPSetType(ksp, KSPPREONLY);
				PCSetType(pc, "lu");
				PCFactorSetMatSolverPackage(pc, "mumps");
				KSPSetInitialGuessNonzero(ksp, PETSC_FALSE);
				KSPSetTolerances(ksp, 0, 0, PETSC_DEFAULT, 1);
			}
			
			if(this->verbose()) {
				SNESMonitorSet(
				snes,
				[](SNES, PetscInt iter, PetscReal res, void*) -> PetscErrorCode {
					PrintInfo::print_iter_status({static_cast<Scalar>(iter), res}); 
					return 0;
				},
				nullptr,
				nullptr);
			}

			// printf("%g %g %g %d\n", this->atol(), this->rtol(), this->stol(), this->max_it());
			SNESSetTolerances(snes, this->atol(), this->rtol(), this->stol(), this->max_it(), 1000);

			SNESLineSearch linesearch;
			SNESGetLineSearch(snes, &linesearch);
			SNESLineSearchSetFromOptions(linesearch);
			SNESLineSearchSetType(linesearch, line_search_type_);

#if !UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)  
			SNESLineSearchSetTolerances(linesearch, PETSC_DEFAULT, PETSC_DEFAULT, this->rtol(), this->atol(), PETSC_DEFAULT, 20);
// 			if(std::string(SNESLINESEARCHBASIC) == line_search_type_) {
// 				SNESLineSearchSetComputeNorms(linesearch, PETSC_FALSE);
// 			}
#endif
			
			if(this->verbose()) {
				this->init_solver("utopia/petsc SemismoothNewton",  {" it.", "|| Au - b||"});
			}

			SNESSolve(snes, nullptr, raw_type(x));

			if(this->verbose()) {
				SNESConvergedReason  reason;
				PetscInt            its; 
				SNESGetConvergedReason(snes, &reason);
				SNESGetIterationNumber(snes, &its);

				this->exit_solver(its, reason); 
			}

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

		void init_snes(SNES &snes)
		{ }

		void set_line_search_type(SNESLineSearchType line_search_type)
		{
			line_search_type_ = line_search_type;
		}
		
	private:	
		std::shared_ptr <Solver>        linear_solver_;
		BoxConstraints                  constraints_;
		SNESLineSearchType line_search_type_;

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
			// PetscReal mag_f = 0.;
			// VecNorm(f, NORM_2, &mag_f);
			// std::cout << "mag_f: " << mag_f << std::endl;
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

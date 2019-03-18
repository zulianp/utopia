#ifndef UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP
#define UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP

#include "utopia_petsc_KSPSolver.hpp"

#include "utopia_Preconditioner.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_SolutionStatus.hpp"
#include "utopia_petsc_KSPTypes.hpp"
#include "utopia_make_unique.hpp"


#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petsc/private/kspimpl.h>

namespace utopia {

	class KSPLog {
	public:
		Vec x_k_1;
		Vec x_k_2;

		KSPLog()
		: x_k_1(nullptr), x_k_2(nullptr)
		{}

		inline bool initialized() const
		{
			return x_k_1 != nullptr;
		}

		void init_from(Vec x)
		{
			destroy();

			VecDuplicate(x, &x_k_1);
			VecDuplicate(x, &x_k_2);
		}

		~KSPLog()
		{
			destroy();
		}

		inline void destroy()
		{
			if(x_k_1) {
				VecDestroy(&x_k_1);
				x_k_1 = nullptr;
			}

			if(x_k_2) {
				VecDestroy(&x_k_2);
				x_k_2 = nullptr;
			}
		}
	};

	template<typename Matrix, typename Vector>
	class KSPSolver<Matrix, Vector, PETSC>::Impl{
	public:
	        Impl(MPI_Comm comm)
	        : ksp_(nullptr), owner_(true)
	        {
	            init(comm);

	            ksp_type(KSPTypes::instance().ksp(0));
	            pc_type(KSPTypes::instance().pc(0));
	            solver_package(KSPTypes::instance().package(0));
	        }

	        Impl(KSP &ksp, const bool owner = false)
	        : ksp_(ksp), owner_(owner)
	        {
	        	set_command_line_parameters(); 
	        }

	        /* @brief      Sets the choice of direct solver.
	         *             Please note, in petsc, direct solver is used as preconditioner alone, with proper settings.
	         *
	         * @param[in]  PCType  The type of direct solver.
	         */
	        inline void pc_type(const std::string &pc_type)
	        {
	            if(KSPTypes::instance().is_pc_valid(pc_type)) {
	                PC pc;
	                KSPGetPC(ksp_, &pc);
	                PCSetType(pc, pc_type.c_str());
	            }
	        }

	        /**
	         * @brief      Sets KSP type
	         */
	        inline void ksp_type(const std::string & ksp_type)
	        {
	            if(KSPTypes::instance().is_ksp_valid(ksp_type)) {
	                KSPSetType(ksp_, ksp_type.c_str());
	            }
	        }

	        /**
	         * @brief      Sets solver package for choice of direct solver.
	         *
	         * @param[in]  SolverPackage  The solver package.
	         */
	        inline void solver_package(const std::string &package)
	        {
	            if(KSPTypes::instance().is_solver_package_valid(package)) {

	            	PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            	PC pc;
	            	ierr = KSPGetPC(ksp_, &pc);  assert(ierr == 0);

	#if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
	                ierr = PCFactorSetMatSolverPackage(pc, package.c_str()); assert(ierr == 0);
	#else
	                ierr = PCFactorSetMatSolverType(pc, package.c_str());    assert(ierr == 0);
	#endif
	            } else {
	            	std::cerr << "Invalid solver package \"" << package << "\"" << std::endl;
	            }
	        }


	        /**
	         * @brief      Returns type of direct solver.
	         */
	        inline PCType pc_type() const
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            PC pc;
	            PCType ret;
	            ierr = KSPGetPC(ksp_, &pc); assert(ierr == 0);
	            ierr = PCGetType(pc, &ret); assert(ierr == 0);
	            return ret;
	        }

	        inline bool has_shell_pc() const
	        {
	            return std::string(PCSHELL) == pc_type();
	        }

	        /**
	         * @brief      Returns ksp package type
	         */
	        inline KSPType ksp_type() const
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            KSPType ret;
	            ierr = KSPGetType(ksp_, &ret); assert(ierr == 0);
	            return ret;
	        }

	        /**
	         * @brief      Returns type of solver package.
	         */
	        inline std::string solver_package() const
	        {
	
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            
	            PC pc;
	            ierr = KSPGetPC(ksp_, &pc);                     assert(ierr == 0);

#if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
	            const MatSolverPackage stype;
	            ierr = PCFactorGetMatSolverPackage(pc, &stype); assert(ierr == 0);
	            if(stype) return stype; else return "";
#else
	            MatSolverType stype;
	            ierr = PCFactorGetMatSolverType(pc, &stype); assert(ierr == 0);
	            if(stype) return stype; else return "";
#endif
	            
	            // static const char * ret = " ";
	            // return ret;
	
	        }

	        inline  bool is_null() const
	        {
	            return ksp_ == nullptr;
	        }

	        ~Impl()
	        {
	            destroy();
	        }

	        inline void destroy()
	        {
	            if(ksp_) {
	                if(owner_) {
	                    KSPDestroy(&ksp_);
	                }

	                ksp_ = nullptr;
	            }
	        }

	        inline MPI_Comm communicator() const {
	            MPI_Comm comm = PetscObjectComm((PetscObject)ksp_);
	            assert(comm != MPI_COMM_NULL);
	            return comm;
	        }

	        void describe() const
	        {
	            PetscBool bool_value;
	            KSPGetInitialGuessNonzero(ksp_, &bool_value);
	            std::cout << "-------------------------------------------" << std::endl;
	            std::cout << "KSPGetInitialGuessNonzero: " << bool_value << std::endl;
	            std::cout << "PCType:                    " << pc_type() << std::endl;
	            std::cout << "KSPType:                   " << ksp_type() << std::endl;
	            std::cout << "-------------------------------------------" << std::endl;

	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            ierr = KSPView(ksp_, PETSC_VIEWER_STDOUT_(communicator())); assert(ierr == 0);
	        }

	        inline void init(const MPI_Comm comm)
	        {
	            destroy();
	            KSPCreate(comm, &ksp_);
	            KSPSetFromOptions(ksp_);
	            KSPSetComputeSingularValues(ksp_, PETSC_FALSE);
	            // KSPSetNormType(ksp_, KSP_NORM_UNPRECONDITIONED);
	        }

	        inline void set_initial_guess_non_zero(const bool val)
	        {
	            if(val) {
	                KSPSetInitialGuessNonzero(ksp_, PETSC_TRUE);
	            } else {
	                KSPSetInitialGuessNonzero(ksp_, PETSC_FALSE);
	            }
	        }

	        void pc_copy_settings(PC &other_pc) const
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            PC this_pc;
	            PCType type;

	            ierr = KSPGetPC(ksp_, &this_pc);    assert(ierr == 0);

	            ierr = PCGetType(this_pc, &type);   assert(ierr == 0);
	            ierr = PCSetType(other_pc, type);   assert(ierr == 0);

#if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
	            const MatSolverPackage stype;
	            ierr = PCFactorGetMatSolverPackage(this_pc, &stype); assert(ierr == 0);
	            ierr = PCFactorSetMatSolverPackage(other_pc, stype); assert(ierr == 0);
#else	
	            MatSolverType stype;
	            ierr =  PCFactorGetMatSolverType(this_pc, &stype); assert(ierr == 0);
	            ierr =  PCFactorSetMatSolverType(other_pc, stype); assert(ierr == 0);
#endif
	        }

	        void copy_settings_from(const Impl &other)
	        {
	            other.copy_settings_to(ksp_);
	        }

	        void copy_settings_to(KSP &other_ksp) const
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            PetscBool bool_value;
	            KSPType type;
	            PetscReal rtol;
	            PetscReal abstol;
	            PetscReal dtol;
	            PetscInt maxits;
	            PC other_pc;

	            KSPSetFromOptions(other_ksp);

	            ierr = KSPGetComputeSingularValues(ksp_, &bool_value);     assert(ierr == 0);
	            ierr = KSPSetComputeSingularValues(other_ksp, bool_value); assert(ierr == 0);

	            ierr = KSPGetType(ksp_, &type);                            assert(ierr == 0);
	            ierr = KSPSetType(other_ksp, type);                        assert(ierr == 0);

	            ierr = KSPGetInitialGuessNonzero(ksp_, &bool_value);       assert(ierr == 0);
	            ierr = KSPSetInitialGuessNonzero(other_ksp, bool_value);   assert(ierr == 0);

	            ierr = KSPGetTolerances(ksp_, &rtol, &abstol, &dtol, &maxits);  assert(ierr == 0);
	            ierr = KSPSetTolerances(other_ksp, rtol, abstol, dtol, maxits); assert(ierr == 0);

	            ierr = KSPGetPC(other_ksp, &other_pc);

	            pc_copy_settings(other_pc);
	        }

	        void attach_shell_preconditioner(
	            PetscErrorCode (*apply)(PC, Vec, Vec),
	            void *ctx,
	            PetscErrorCode (*setup)(PC),
	            PetscErrorCode (*destroy)(PC)
	            )
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            PC pc;

	            ierr = KSPGetPC(ksp_, &pc);                         assert(ierr == 0);
	            ierr = PCSetType(pc, PCSHELL);                      assert(ierr == 0);
	            ierr = PCShellSetApply(pc, apply);                  assert(ierr == 0);

	            ierr = PCShellSetName(pc, "Utopia Preconditioner"); assert(ierr == 0);

	            if(setup) {
	                ierr = PCShellSetSetUp(pc, setup);              assert(ierr == 0);
	            }

	            if(destroy) {
	                ierr = PCShellSetSetUp(pc, destroy);            assert(ierr == 0);
	            }

	            if(ctx) {
	                ierr = PCShellSetContext(pc, ctx);              assert(ierr == 0);
	            }

	            //usefull methods
	            //PCGetOperatorsSet(PC pc,PetscBool  *mat,PetscBool  *pmat)
	            //PetscErrorCode  PCGetOperators(PC pc,Mat *Amat,Mat *Pmat)
	            //PetscErrorCode  PCShellSetApplyRichardson(PC pc,PetscErrorCode (*apply)(PC,Vec,Vec,Vec,PetscReal,PetscReal,PetscReal,PetscInt,PetscBool,PetscInt*,PCRichardsonConvergedReason*))
	            //PCGetReusePreconditioner(PC pc,PetscBool *flag)
	        }

	        void update(const Matrix &mat)
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);

	            ierr = KSPSetOperators(ksp_, raw_type(mat), raw_type(mat)); assert(ierr == 0);
	            ierr = KSPSetUp(ksp_);                                      assert(ierr == 0);
	        }

	        void update(const Matrix &mat, const Matrix &prec)
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);

	            ierr = KSPSetOperators(ksp_, raw_type(mat), raw_type(prec)); assert(ierr == 0);
	            ierr = KSPSetUp(ksp_);                                       assert(ierr == 0);
	        }

	        bool smooth(const SizeType sweeps,
	                    const Vector &rhs,
	                    Vector &x)
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            KSPNormType normtype;
	            PetscReal rtol;
	            PetscReal abstol;
	            PetscReal dtol;
	            PetscInt maxits;
	            void *ctx;

	            //back-up settings
	            ierr = KSPGetNormType(ksp_, &normtype);                                                    assert(ierr == 0);
	            ierr = KSPGetTolerances(ksp_, &rtol, &abstol, &dtol, &maxits);                             assert(ierr == 0);

	            //set up smoothing things
	            ierr = KSPSetTolerances(ksp_, 0., 0., PETSC_DEFAULT, sweeps);                              assert(ierr == 0);
	            ierr = KSPSetNormType(ksp_, KSP_NORM_NONE);                                                assert(ierr == 0);
	            ierr = KSPSetConvergenceTest(ksp_, KSPConvergedSkip, NULL, NULL);                          assert(ierr == 0);

	            //perform smoothing
	            ierr = KSPSolve(ksp_, raw_type(rhs), raw_type(x));                                         assert(ierr == 0);

	            //reset-solver-stuff
	            ierr = KSPSetNormType(ksp_, normtype);                                                     assert(ierr == 0);
	            ierr = KSPSetTolerances(ksp_, rtol, abstol, dtol, maxits);                                 assert(ierr == 0);
	            ierr = KSPConvergedDefaultCreate(&ctx);                                                    assert(ierr == 0);
	            ierr = KSPSetConvergenceTest(ksp_,  KSPConvergedDefault, ctx, KSPConvergedDefaultDestroy); assert(ierr == 0);
	            return true;
	        }

	        void set_tolerances(const PetscReal rtol,
	                            const PetscReal atol,
	                            const PetscReal dtol,
	                            const PetscInt max_it)
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            ierr = KSPSetTolerances(ksp_, rtol, atol, dtol, max_it); assert(ierr == 0);
	        }

	        void number_of_subdomains(const PetscInt &n_blocks)
	        {
	        	PetscErrorCode ierr;
	        	PC pc;
	            PCType pc_type;

	            ierr = KSPGetPC(ksp_, &pc);    assert(ierr == 0);
	            ierr = PCGetType(pc, &pc_type);   assert(ierr == 0);

	            PetscBool flg_bjacobi=PETSC_FALSE, flg_asm=PETSC_FALSE; 
	            PetscObjectTypeCompare((PetscObject)pc,PCBJACOBI,&flg_bjacobi);

	            if(!flg_bjacobi){
	            	PetscObjectTypeCompare((PetscObject)pc,PCASM,&flg_asm);
	            }

	            if(flg_bjacobi)
	            {
	            	PCBJacobiSetTotalBlocks(pc, n_blocks, nullptr); 
	            }
	            else if(flg_asm)
	            {
	            	PCGASMSetTotalSubdomains(pc, n_blocks); 
	            }
	            else
	            {
	            	utopia_error("Number of subdomain can be set only for PCBJACOBI or PCASM."); 
	            }
	        }

	        void overlap(const PetscInt &n_overlap)
	        {
	        	PetscErrorCode ierr;
	        	PC pc;
	            PCType pc_type;

	            ierr = KSPGetPC(ksp_, &pc);    assert(ierr == 0);
	            ierr = PCGetType(pc, &pc_type);   assert(ierr == 0);

	            PetscBool flg_asm=PETSC_FALSE; 
	            PetscObjectTypeCompare((PetscObject)pc,PCASM,&flg_asm);
	            
				if(flg_asm)
	            {
	            	PCGASMSetOverlap(pc, n_overlap); 
	            }
	            else
	            {
	            	utopia_error("Overlap can be set only for PCASM"); 
	            	return; 
	            }	
	        }


	        void sub_ksp_pc_type(const std::string sub_ksp_type, const std::string sub_pc_type)
	        {
	        	if(!ksp_->setupstage)
	        	{
	        		utopia_error("sub_ksp_pc_type can be only called after update(). "); 
	        		return; 
	        	}

	        	PetscErrorCode ierr;
	        	PC pc;
	            PCType pc_type;

	            ierr = KSPGetPC(ksp_, &pc);    assert(ierr == 0);
	            ierr = PCGetType(pc, &pc_type);   assert(ierr == 0);

	            PetscBool flg_bjacobi=PETSC_FALSE, flg_asm=PETSC_FALSE; 
	            PetscObjectTypeCompare((PetscObject)pc,PCBJACOBI,&flg_bjacobi);


	            PetscInt first, nlocal; 
	            KSP *subksp; 

	            if(!flg_bjacobi){
	            	PetscObjectTypeCompare((PetscObject)pc,PCASM,&flg_asm);
	            }

	            if(flg_bjacobi)
	            {
	            	PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
	            }
	            else if(flg_asm)
	            {
	            	PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
	            }
	            else
	            {
	            	utopia_error("Sub_ksp_pc_type can be set only for PCBJACOBI and PCASM"); 
	            	return; 
	            }

				for (auto i=0; i<nlocal; i++) 
			  	{
			  		KSPSetType(subksp[i], sub_ksp_type.c_str());
			  		PC subpc; 
			  	  	KSPGetPC(subksp[i],&subpc);
			  	    PCSetType(subpc, sub_pc_type.c_str());
			  	}
	        }


	        void sub_ksp_type(const std::string sub_ksp_type)
	        {
	        	if(!ksp_->setupstage)
	        	{
	        		utopia_error("sub_ksp_type can be only called after update(). "); 
	        		return; 
	        	}

	        	PetscErrorCode ierr;
	        	PC pc;
	            PCType pc_type;

	            ierr = KSPGetPC(ksp_, &pc);    assert(ierr == 0);
	            ierr = PCGetType(pc, &pc_type);   assert(ierr == 0);

	            PetscBool flg_bjacobi=PETSC_FALSE, flg_asm=PETSC_FALSE; 
	            PetscObjectTypeCompare((PetscObject)pc,PCBJACOBI,&flg_bjacobi);


	            PetscInt first, nlocal; 
	            KSP *subksp; 

	            if(!flg_bjacobi){
	            	PetscObjectTypeCompare((PetscObject)pc,PCASM,&flg_asm);
	            }

	            if(flg_bjacobi)
	            {
	            	PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
	            }
	            else if(flg_asm)
	            {
	            	PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
	            }
	            else
	            {
	            	utopia_error("Sub_ksp_pc_type can be set only for PCBJACOBI and PCASM"); 
	            	return; 
	            }

            	for (auto i=0; i<nlocal; i++) 
			  	{
			  		KSPSetType(subksp[i], sub_ksp_type.c_str());
			  	}	            
	        }

	        void sub_pc_type(const std::string sub_pc_type)
	        {
	        	if(!ksp_->setupstage)
	        	{
	        		utopia_error("sub_ksp_pc_type can be only called after update(). "); 
	        		return; 
	        	}

	        	PetscErrorCode ierr;
	        	PC pc;
	            PCType pc_type;

	            ierr = KSPGetPC(ksp_, &pc);    assert(ierr == 0);
	            ierr = PCGetType(pc, &pc_type);   assert(ierr == 0);

	            PetscBool flg_bjacobi=PETSC_FALSE, flg_asm=PETSC_FALSE; 
	            PetscObjectTypeCompare((PetscObject)pc,PCBJACOBI,&flg_bjacobi);


	            PetscInt first, nlocal; 
	            KSP *subksp; 

	            if(!flg_bjacobi){
	            	PetscObjectTypeCompare((PetscObject)pc,PCASM,&flg_asm);
	            }

	            if(flg_bjacobi)
	            {
	            	PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
	            }
	            else if(flg_asm)
	            {
	            	PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
	            }
	            else
	            {
	            	utopia_error("Sub_ksp_pc_type can be set only for PCBJACOBI and PCASM"); 
	            	return; 
	            }

				for (auto i=0; i<nlocal; i++) 
			  	{
			  		PC subpc; 
			  	  	KSPGetPC(subksp[i],&subpc);
			  	    PCSetType(subpc, sub_pc_type.c_str());
			  	}
	        }	        


	        void sub_solver_package(const std::string sub_pc_type)
	        {
	        	if(!ksp_->setupstage)
	        	{
	        		utopia_error("sub_ksp_pc_type can be only called after update(). "); 
	        		return; 
	        	}

	        	PetscErrorCode ierr;
	        	PC pc;
	            PCType pc_type;

	            ierr = KSPGetPC(ksp_, &pc);    assert(ierr == 0);
	            ierr = PCGetType(pc, &pc_type);   assert(ierr == 0);

	            PetscBool flg_bjacobi=PETSC_FALSE, flg_asm=PETSC_FALSE; 
	            PetscObjectTypeCompare((PetscObject)pc,PCBJACOBI,&flg_bjacobi);


	            PetscInt first, nlocal; 
	            KSP *subksp; 

	            if(!flg_bjacobi){
	            	PetscObjectTypeCompare((PetscObject)pc,PCASM,&flg_asm);
	            }

	            if(flg_bjacobi)
	            {
	            	PCBJacobiGetSubKSP(pc, &nlocal, &first, &subksp);
	            }
	            else if(flg_asm)
	            {
	            	PCASMGetSubKSP(pc, &nlocal, &first, &subksp);
	            }
	            else
	            {
	            	utopia_error("Sub_ksp_pc_type can be set only for PCBJACOBI and PCASM"); 
	            	return; 
	            }

				for (auto i=0; i<nlocal; i++) 
			  	{
			  		PC subpc; 
			  	  	KSPGetPC(subksp[i],&subpc);
			  	    #if UTOPIA_PETSC_VERSION_LESS_THAN(3,9,0)
	                	PCFactorSetMatSolverPackage(subpc, sub_pc_type.c_str()); 
					#else
	                	PCFactorSetMatSolverType(subpc, sub_pc_type.c_str());  
					#endif
			  	}
	        }	        

	        void set_monitor(
	            PetscErrorCode (*monitor)(KSP,PetscInt,PetscReal,void *),
	            void *mctx,
	            PetscErrorCode (*monitordestroy)(void**)
	        )
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            ierr = KSPMonitorSet(ksp_, monitor, mctx, monitordestroy); assert(ierr == 0);
	        }

	        bool apply(const Vector &b, Vector &x)
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            auto ls = local_size(b).get(0);
	            auto gs = size(b).get(0);

	            if(gs != size(x).get(0))
	            {
	                x = local_zeros(ls);
	            }

	            assert(ls == local_size(x).get(0));

	            KSPConvergedReason  reason;
	            ierr = KSPSolve(ksp_, raw_type(b), raw_type(x)); assert(ierr == 0);
	            ierr = KSPGetConvergedReason(ksp_, &reason);     assert(ierr == 0);

	            if(reason < 0) {

	                utopia_warning(
	                	"ksp apply returned " + std::to_string(reason) + " = " + converged_str(reason) + 
	                	" ksp_type=" + ksp_type() + " pc_type=" + pc_type() + " solver_package=" + solver_package());
	            }

	            return reason >= 0;
	        }

	        void solution_status(SolutionStatus &status)
	        {
	            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
	            PetscInt its;
	            KSPConvergedReason reason;
	            PetscReal rnorm;

	            ierr = KSPGetIterationNumber(ksp_, &its);    assert(ierr == 0);
	            ierr = KSPGetConvergedReason(ksp_, &reason); assert(ierr == 0);
	            ierr = KSPGetResidualNorm(ksp_, &rnorm);     assert(ierr == 0);

	            status.iterates = its;
	            status.reason = reason;
	            status.gradient_norm = rnorm;
	        }

	        void set_command_line_parameters()
	        {
	            // petsc command line options
	            char           name_[1024];
	            PetscBool      flg;

	#if UTOPIA_PETSC_VERSION_LESS_THAN(3,7,0)
	            PetscOptionsGetString(nullptr, "-ksp_type", name_, 1024, &flg);
	            if(flg) {
	                ksp_type(name_);
	            }

	            PetscOptionsGetString(nullptr, "-pc_type", name_, 1024, &flg);
	            if(flg) {
	                pc_type(name_);
	            }

	            PetscOptionsGetString(nullptr, "-pc_factor_mat_solver_package", name_, 1024, &flg);
	            if(flg) {
	                solver_package(name_);
	            }
	        
	#else
	            PetscOptionsGetString(nullptr, nullptr, "-ksp_type", name_, 1024, &flg);
	            if(flg) {
	                ksp_type(name_);
	            }

	            PetscOptionsGetString(nullptr, nullptr, "-pc_type", name_, 1024, &flg);
	            if(flg) {
	                pc_type(name_);
	            }

	            PetscOptionsGetString(nullptr, nullptr, "-pc_factor_mat_solver_package", name_, 1024, &flg);
	            if(flg) {
	                solver_package(name_);
	            }
	            
	          	PetscOptionsGetString(nullptr, nullptr, "-pc_factor_mat_solver_type", name_, 1024, &flg);
	          	if(flg) {
	             	 solver_package(name_);
	          	}
	#endif
	        }

	        KSP &implementation()
	        {
	            return ksp_;
	        }

	        const KSP &implementation() const
	        {
	            return ksp_;
	        }

	    private:
	        KSP            ksp_;
	        bool owner_;
	};


	///////////////////////////////////////////////////////////////////////////////////////////////////

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::KSPSolver()
	: ksp_(utopia::make_unique<Impl>(PETSC_COMM_WORLD))
	{
		ksp_type("bcgs");
		pc_type("jacobi");
		ksp_->set_initial_guess_non_zero(true);
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::KSPSolver(std::unique_ptr<typename utopia::KSPSolver<Matrix, Vector, PETSC>::Impl> &&w)
	: ksp_(std::move(w))
	{}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::wrap(KSP &ksp)
	{
		ksp_ = utopia::make_unique<Impl>(ksp, false);
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::~KSPSolver() {}


	template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_initial_guess_non_zero(const bool val)
	{
		ksp().set_initial_guess_non_zero(val);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::pc_type(const std::string &pc_type)
	{
		ksp_->pc_type(pc_type);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::ksp_type(const std::string & ksp_type)
	{
		ksp_->ksp_type(ksp_type);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::solver_package(const std::string &package)
	{
		ksp_->solver_package(package);
	}


    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::number_of_subdomains(const SizeType & n)
	{
		return this->ksp_->number_of_subdomains(n);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::overlap(const SizeType & n)
	{
		return this->ksp_->overlap(n);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::sub_ksp_pc_type(const std::string ksp_type, const std::string pc_type)
	{
		return this->ksp_->sub_ksp_pc_type(ksp_type, pc_type); 
	}
        
    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::sub_ksp_type(const std::string type)
	{
		return this->ksp_->sub_ksp_type(type); 
	} 

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::sub_pc_type(const std::string type)
	{
		return this->ksp_->sub_pc_type(type); 
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::sub_solver_package(const std::string type)
	{
		return this->ksp_->sub_solver_package(type); 
	} 

    template<typename Matrix, typename Vector>
	std::string KSPSolver<Matrix, Vector, PETSC>::pc_type() const { return this->ksp_->pc_type();}

    template<typename Matrix, typename Vector>
	std::string KSPSolver<Matrix, Vector, PETSC>::ksp_type() const { return this->ksp_->ksp_type();}

    template<typename Matrix, typename Vector>
	std::string KSPSolver<Matrix, Vector, PETSC>::solver_package() const { return this->ksp_->solver_package();}

    template<typename Matrix, typename Vector>
	bool KSPSolver<Matrix, Vector, PETSC>::apply(const Vector &b, Vector &x)
	{
		ksp_->set_tolerances(this->rtol(), this->atol(), PETSC_DEFAULT, this->max_it());

        // is this proper place to do so???
        // this->set_ksp_options(ksp_->implementation());
		return ksp_->apply(b, x);
	}

    template<typename Matrix, typename Vector>
	bool KSPSolver<Matrix, Vector, PETSC>::smooth(const Vector &rhs, Vector &x)
	{
		return ksp_->smooth(this->sweeps(), rhs, x);
	}

    template<typename Matrix, typename Vector>
	bool KSPSolver<Matrix, Vector, PETSC>::must_compute_cond_number() const
	{
		static const std::vector<std::string> types = {"bcgs", "cg", "gmres"};
		return (this->verbose() && in_array(ksp_->ksp_type(), types));
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_ksp_options(KSP &ksp)
	{
		ksp_->copy_settings_to(ksp);
		set_monitor_options(ksp);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::attach_preconditioner(KSP &ksp) const {
		Impl w(ksp, false);
		auto delegate_ptr = std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());

		if(delegate_ptr) {
			if(ksp_->has_shell_pc()) {
				m_utopia_warning_once("set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
			}
		} else if(this->get_preconditioner()) {
			auto shell_ptr = this->get_preconditioner().get();
			w.attach_shell_preconditioner(UtopiaPCApplyShell,
				shell_ptr,
				nullptr,
				nullptr
				);
		}
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
	{
		PreconditionedSolver::set_preconditioner(precond);

		auto delegate_ptr = std::dynamic_pointer_cast<DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());

		if(delegate_ptr) {
			if(ksp_->has_shell_pc()) {
				m_utopia_warning_once("set_preconditioner sets jacobi if a delegate precond has been set and type is matshell");
				ksp_->pc_type("jacobi");
			}

		} else if(this->get_preconditioner()) {
			auto shell_ptr = this->get_preconditioner().get();
			ksp_->attach_shell_preconditioner(UtopiaPCApplyShell,
				shell_ptr,
				nullptr,
				nullptr
				);
		}
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::set_monitor_options(KSP &ksp) const
	{
		PetscErrorCode ierr; UTOPIA_UNUSED(ierr);
		if(this->verbose()) {
			ierr = KSPMonitorSet(ksp,
				MyKSPMonitor,
				new KSPLog(),
				[](void **arg) -> PetscErrorCode {
					auto ksp_log = (KSPLog **)(arg);
					delete *ksp_log;
					*ksp_log = nullptr;
					return 0;
				});  assert(ierr == 0);
		}

		if(must_compute_cond_number()) {
			ierr = KSPSetComputeSingularValues(ksp, PETSC_TRUE); assert(ierr == 0);
		}
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::handle_reset(const Matrix &op)
	{
		bool must_reset = true;
        // bool must_reset = false;

        // if(this->has_operator()) {
        //     auto s = size(*this->get_operator());
        //     auto other_s = size(op);
        //     if(s != other_s) {
        //     must_reset = true;
        //     }
        // }

        // if(op.implementation().communicator() != ksp_->communicator())
        // {
        //     must_reset = true;
        // }

		if(must_reset) {
			auto temp_ksp = utopia::make_unique<Impl>(op.implementation().communicator());
			temp_ksp->copy_settings_from(*ksp_);
			ksp_ = std::move(temp_ksp);

			if(this->get_preconditioner()) {
				set_preconditioner(this->get_preconditioner());
			}
		}
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op, const std::shared_ptr<const Matrix> &prec)
	{
		handle_reset(*op);
		set_monitor_options(ksp_->implementation());


		PreconditionedSolver::update(op, prec);
		ksp_->update(*op, *prec);
	}

    template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::update(const std::shared_ptr<const Matrix> &op)
	{
		handle_reset(*op);

		PreconditionedSolver::update(op);
		set_monitor_options(ksp_->implementation());

		bool skip_set_operators = false;
		if(this->get_preconditioner()) {
			auto mat_prec = std::dynamic_pointer_cast< DelegatePreconditioner<Matrix, Vector>>(this->get_preconditioner());
			if(mat_prec) {
				ksp_->update(*op, *mat_prec->get_matrix());
				skip_set_operators = true;
			}
		}

		if(!skip_set_operators) {
			ksp_->update(*op);
		}
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC> & KSPSolver<Matrix, Vector, PETSC>::operator=(const KSPSolver<Matrix, Vector, PETSC> &other)
	{
		if(this == &other) return *this;
		PreconditionedSolver::operator=(other);
		Smoother::operator=(other);

		ksp_ = utopia::make_unique<Impl>(other.ksp_->communicator());
		ksp_->copy_settings_from(*other.ksp_);
		return *this;
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC> & KSPSolver<Matrix, Vector, PETSC>::operator=(KSPSolver<Matrix, Vector, PETSC> &&other)
	{
		if(this == &other) return *this;
		PreconditionedSolver::operator=(std::move(other));
		Smoother::operator=(std::move(other));
		ksp_ = std::move(other.ksp_);
		return *this;
	}


    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::KSPSolver(const KSPSolver<Matrix, Vector, PETSC> &other):
	PreconditionedSolver(other),
	Smoother(other),
	ksp_(utopia::make_unique<Impl>(other.ksp_->communicator()))
	{
		ksp_->copy_settings_from(*other.ksp_);
	}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC>::KSPSolver(KSPSolver<Matrix, Vector, PETSC> &&other)
	: PreconditionedSolver(std::move(other)),
	Smoother(std::move(other)),
	ksp_(std::move(other.ksp_))
	{}

    template<typename Matrix, typename Vector>
	KSPSolver<Matrix, Vector, PETSC> * KSPSolver<Matrix, Vector, PETSC>::clone() const
	{
		return new KSPSolver(*this);
	}

    template<typename Matrix, typename Vector>
	typename utopia::KSPSolver<Matrix, Vector, PETSC>::Impl &KSPSolver<Matrix, Vector, PETSC>::ksp()
	{
		return *ksp_;
	}

    template<typename Matrix, typename Vector>
	const typename utopia::KSPSolver<Matrix, Vector, PETSC>::Impl &KSPSolver<Matrix, Vector, PETSC>::ksp() const
	{
		return *ksp_;
	}

	template<typename Matrix, typename Vector>
	KSP &KSPSolver<Matrix, Vector, PETSC>::implementation()
	{
		return ksp().implementation();
	}

	template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::describe(std::ostream &os) const
	{
		os << "KSPSolver: \n";
		os << "ksp_type:       " << ksp_type() 		 << "\n";
		os << "pc_type:        " << pc_type() 		 << "\n";
		os << "solver_package: " << solver_package() << "\n";
	}

	template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::read(Input &is)
	{
		PreconditionedSolver::read(is);

		//TODO add other options
		auto new_ksp_type = this->ksp_type();
		auto new_pc_type = this->pc_type();
		auto new_solver_package = this->solver_package();


		is.get("ksp_type", new_ksp_type);
		is.get("pc_type", new_pc_type);
		is.get("solver_package", new_solver_package);

		this->ksp_type(new_ksp_type);
		this->pc_type(new_pc_type);
		this->solver_package(new_solver_package);
	}

	template<typename Matrix, typename Vector>
	void KSPSolver<Matrix, Vector, PETSC>::print_usage(std::ostream &os) const
	{
		PreconditionedSolver::print_usage(os);
		this->print_param_usage(os, "ksp_type", "string", "Type of KSP solver.", "bcgs"); 
		this->print_param_usage(os, "pc_type", "string", "Type of preconditoner.", "jacobi"); 
		this->print_param_usage(os, "solver_package", "string", "Type of solver package.", " "); 
	}
}

#endif //UTOPIA_PETSC_KSP_SOLVER_IMPL_HPP

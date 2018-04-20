/*
* @Author: Alena Kopanicakova
* @Date:   2016-09-01
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-06-09
*/
#ifndef UTOPIA_PETSC_KSP_H
#define UTOPIA_PETSC_KSP_H

#include "utopia_Preconditioner.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_Smoother.hpp"

#include "utopia_Core.hpp"

#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>

namespace utopia {


    PetscErrorCode UtopiaPCApplyShell(PC pc, Vec x, Vec y);
    PetscErrorCode MyKSPMonitor(KSP,PetscInt,PetscReal,void*);

    typedef struct
    {
        Vec x_k_1;       
        Vec x_k_2;       

        PetscBool compute_cond_number; 
    }
    UTOPIA_TRACE;

    
    /**@ingroup     Linear 
     * @brief       Class provides interface to Petsc KSP solvers \n
     *              For setting up basic parameters, one can use classic Petsc runtime options, e.g. 
     *              To see all possibilities, please refer to: 
     *                                                   * <a href="http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html">preconditioner types</a> 
     *                                                   * <a href="http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType">solver types</a>
     *              
     *              Setting own/utopia preconditioner, can be done as following: 
     *              \snippet tests/utopia_SolverTest.cpp PetscKSPSolver solve example1
     *              Detailed information about preconditioners, can be found in  \ref precondotioners.
     */
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
    class KSPSolver {};


    template<typename Matrix, typename Vector>
    class KSPSolver<Matrix, Vector, PETSC> : virtual public PreconditionedSolver<Matrix, Vector>, public Smoother<Matrix, Vector>
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver; 
        typedef utopia::PreconditionedSolver<Matrix, Vector> PreconditionedSolver;

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "only works with petsc types");



        KSPSolver(  const Parameters params = Parameters(), 
                    const std::vector<std::string> ksp_types    = {"bcgs", "cg", "groppcg", "pipecg", "pipecgrr", "fcg", "pipefcg", "gmres", "pipefgmres",   "fgmres",   "lgmres",   "dgmres",   "pgmres", "tcqmr", "ibcgs",   "fbcgs",   "fbcgsr",   "bcgsl", "cgs", "tfqmr", "cr", "pipecr", "lsqr", "preonly", "qcg", "bicg", "minres", "symmlq", "lcd", "python", "gcr", "pipegcr", "tsirm", "cgls"},
                    const std::vector<std::string> pc_types     = {"jacobi","sor","lu","bjacobi","eisenstat","ilu","icc","asm","gasm","ksp","cholesky","pbjacobi","mat","hypre", "cp","bfbt","lsc","python","pfmg","syspfmg","redistribute","svd","gamg","bicgstabcusp","ainvcusp","bddc"}, 
                    const std::vector<std::string> pc_packages = {" "}):

                KSP_types(ksp_types),
                PC_types(pc_types),
                Solver_packages(pc_packages), 
                compute_cond_number(PETSC_FALSE)

        {
            set_parameters(params); 
            KSP_type_       = KSP_types.at(0); 
            PC_type_        = PC_types.at(0);
            solver_package_ = Solver_packages.at(0);

            //KSPCreate(PETSC_COMM_WORLD, &ksp);  - something is wrong here ... 
        }


        virtual ~KSPSolver()
        { 
            // KSPDestroy(&ksp);
        }


    /**
     * @brief      Sets the parameters.
     *
     * @param[in]  params  The parameters
     */
    virtual void set_parameters(const Parameters params) override
    {
        PreconditionedSolver::set_parameters(params);

        // our param options - useful for passing from passo 
        ksp_type(params.lin_solver_type());                               
        pc_type(params.preconditioner_type());      
        solver_package(params.preconditioner_factor_mat_solver_package());                             
    
        // petsc command line options 
        char           name_[1024]; 
        PetscBool      flg;

        #if UTOPIA_PETSC_VERSION_LESS_THAN(3,7,0)
             PetscOptionsGetString(NULL, "-ksp_type", name_, 1024, &flg);
            if(flg)
                ksp_type(name_);

            PetscOptionsGetString(NULL, "-pc_type", name_, 1024, &flg);
            if(flg)
                pc_type(name_);

            PetscOptionsGetString(NULL, "-pc_factor_mat_solver_package", name_, 1024, &flg);
            if(flg)
                solver_package(name_);
        #else
             PetscOptionsGetString(NULL, NULL, "-ksp_type", name_, 1024, &flg);
            if(flg)
                ksp_type(name_);

            PetscOptionsGetString(NULL, NULL, "-pc_type", name_, 1024, &flg);
            if(flg)
                pc_type(name_);

            PetscOptionsGetString(NULL,NULL, "-pc_factor_mat_solver_package", name_, 1024, &flg);
            if(flg)
                solver_package(name_);
            
        #endif
 
    }




     /* @brief      Sets the choice of direct solver. 
     *             Please note, in petsc, direct solver is used as preconditioner alone, with proper settings. 
     *
     * @param[in]  PCType  The type of direct solver. 
     */
    virtual void pc_type(const std::string & PCType ) 
    { 
        PC_type_ = in_array(PCType, PC_types) ? PCType : PC_types.at(0);
    }; 

    /**
     * @brief      Sets KSP type 
     */
    void ksp_type(const std::string & KSPType ) 
    { 
        KSP_type_ = in_array(KSPType, KSP_types) ? KSPType : KSP_types.at(0);
    }; 

    /**
     * @brief      Sets solver package for choice of direct solver. 
     *
     * @param[in]  SolverPackage  The solver package.
     */
    void solver_package(const std::string & SolverPackage ) 
    { 
        solver_package_ = in_array(SolverPackage, Solver_packages) ? SolverPackage : Solver_packages.at(0); 
    }; 


    /**
     * @brief      Returns type of direct solver. 
     */
    const std::string &pc_type() const { return PC_type_;}; 

    /**
     * @brief      Returns ksp package type 
     */
    const std::string &ksp_type() const { return KSP_type_; }; 


    /**
     * @brief      Returns type of solver package. 
     */
    const std::string &solver_package() const { return solver_package_; }; 



public:

    /**
    * @brief        Solve function for all PETS KSP solvers. 
     *              It is also compatible with our own Utopia preconditioners.
    */
    bool apply(const Vector &b, Vector &x) override
    {
        PetscErrorCode ierr;

        Size ls = local_size(b);
        Size gs = size(b);

        if(empty(x) || gs.get(0) != size(x).get(0)) 
        {
            x = local_zeros(ls.get(0));
        }

        assert(size(b).get(0) == size(x).get(0));
        assert(local_size(b).get(0) == local_size(x).get(0));

        const Matrix &A = *this->get_operator(); 


        PreconditionedSolver::init_solver("UTOPIA::Petsc KSP", {}); 
          
        KSPConvergedReason  reason;
        // PetscReal           r_norm;
        PetscInt            its; 
        MPI_Comm            comm; 
          

        ierr = PetscObjectGetComm((PetscObject)raw_type(A), &comm);
        ierr = KSPCreate(comm, &ksp);

        bool skip_set_operators = false;
        
        if(this->get_preconditioner()) {
            auto mat_prec = dynamic_cast< DelegatePreconditioner<Matrix, Vector> *>(this->get_preconditioner().get());
            if(mat_prec) {
                ierr = KSPSetOperators(ksp, raw_type(A), raw_type(*mat_prec->get_matrix()));
                skip_set_operators = true;
            } 
        }
       
        if(!skip_set_operators) { 
            ierr = KSPSetOperators(ksp, raw_type(A), raw_type(A));
        }

        set_ksp_options(ksp); 

        //  Set the user-defined routine for applying the preconditioner 
        if(!skip_set_operators)
            attach_preconditioner(ksp); 


        VecDuplicate(raw_type(x), &(ut_log.x_k_2));
        VecDuplicate(raw_type(x), &(ut_log.x_k_1));

        ut_log.compute_cond_number = compute_cond_number; 

        ierr = KSPSetUp(ksp);
        ierr = KSPSolve(ksp, raw_type(b), raw_type(x));
          
        ierr = KSPGetConvergedReason(ksp, &reason);
        ierr = KSPGetIterationNumber(ksp, &its);
        
        this->exit_solver(its, reason); 
          
        VecDestroy(&(ut_log.x_k_1)); 
        VecDestroy(&(ut_log.x_k_2)); 

        KSPDestroy(&ksp);

        return true;
    }


    bool smooth(const Matrix &A, const Vector &rhs, Vector &x) override
    {
        KSP solver;
        KSPCreate(A.implementation().communicator(), &solver);
        KSPSetFromOptions(solver); 
        KSPSetType(solver, KSP_type_.c_str());
        KSPSetInitialGuessNonzero(solver, PETSC_TRUE);
        KSPSetTolerances(solver, 0., 0., PETSC_DEFAULT, this->sweeps());

        KSPSetOperators(solver, raw_type(A), raw_type(A));

        if(!this->get_preconditioner()) 
        {
            PC pc; 
            KSPGetPC(solver, &pc);
            PCSetType(pc, PC_type_.c_str());
        }
        
        if(this->verbose()) {
            KSPMonitorSet(
                solver,
                [](KSP, PetscInt iter, PetscReal res, void*) -> PetscErrorCode {
                    PrintInfo::print_iter_status({static_cast<PetscReal>(iter), res}); 
                    return 0;
                },
                nullptr,
                nullptr);
        }
        
        KSPSetUp(solver);
        KSPSolve(solver, raw_type(rhs), raw_type(x));
        KSPDestroy(&solver);
        return true;
    }

    virtual void attach_preconditioner(KSP & ksp)
    {
        PetscErrorCode ierr;
        if(this->get_preconditioner()) 
        {
            PC pc; 
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCSHELL);  
            ierr = PCShellSetApply(pc, UtopiaPCApplyShell);
  
            auto shell_ptr = dynamic_cast<LinearSolver *>(this->get_preconditioner().get());
            ierr = PCShellSetContext(pc, shell_ptr);
            ierr = PCShellSetName(pc,"Utopia Preconditioner");
        }
    }

    /**
     * @brief      Sets the default options for PETSC KSP solver. \n
     *             Default: BiCGstab
     *
     * @param      ksp   The ksp
     */
    virtual void set_ksp_options(KSP & ksp)
    {
        PetscErrorCode ierr;


        // check if our options overwrite this 
        KSPSetFromOptions(ksp); 

        // TODO:: extend to other supported types... 
        compute_cond_number = ( (   KSP_type_ == "bcgs" || KSP_type_ == "cg"  || 
                                    KSP_type_ == "gmres") && (this->verbose() )) ? PETSC_TRUE : PETSC_FALSE; 

        if(this->verbose())
            KSPMonitorSet(ksp, MyKSPMonitor, &ut_log, 0);

        if(compute_cond_number)
            ierr = KSPSetComputeSingularValues(ksp, PETSC_TRUE); 
        
        ierr = KSPSetType(ksp, KSP_type_.c_str());


        // KSP_NORM_PRECONDITIONED  is default in petsc
        //Why doesn't it work with cg?????????????
        if(KSP_type_ != "cg") {
            ierr = KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED); 
        }

        if(!this->get_preconditioner()) 
        {
            PC pc; 
            ierr = KSPGetPC(ksp, &pc);
            ierr = PCSetType(pc, PC_type_.c_str());
        }

        ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
        ierr = KSPSetTolerances(ksp, PreconditionedSolver::rtol(), PreconditionedSolver::atol(), PETSC_DEFAULT,  PreconditionedSolver::max_it());
    }


protected:

    std::string KSP_type_;                                  /*!< Choice of preconditioner types. */  
    const std::vector<std::string> KSP_types;              /*!< Valid options for direct solver types. */  

    std::string  PC_type_;                                  /*!< Choice of preconditioner types. */  
    const std::vector<std::string> PC_types;      

    std::string solver_package_;                            /*!< Choice of direct solver. */     
    const std::vector<std::string> Solver_packages;       /*!< Valid options for Solver packages types. */

    KSP                 ksp;
    UTOPIA_TRACE          ut_log; 
    PetscBool           compute_cond_number;  

    };
    
}



#endif //UTOPIA_PETSC_KSP_H
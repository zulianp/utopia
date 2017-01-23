/*
* @Author: Alena Kopanicakova
* @Date:   2016-09-01
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2017-01-21
*/
#ifndef UTOPIA_PETSC_KSP_H
#define UTOPIA_PETSC_KSP_H

#include "utopia_Core.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_PETSc.hpp"
#include "utopia_Preconditioner.hpp"

#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>

namespace utopia 
{


    PetscErrorCode UtopiaPCApplyShell(PC pc, Vec x, Vec y);
    PetscErrorCode MyKSPMonitor(KSP,PetscInt,PetscReal,void*);

    typedef struct
    {
        Vec x_k_1;       
        Vec x_k_2;       
    }
    UTOPIA_LOG;

    
    /**@ingroup     Linear 
     * @brief       Class provides interface to Petsc KSP solvers \n
     *              For setting up basic parameters, one can use classic PETSc runtime options, e.g. 
     *              To see all possibilities, please refer to: 
     *                                                   * <a href="http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/PC/PCType.html">preconditioner types</a> 
     *                                                   * <a href="http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/KSP/KSPType.html#KSPType">solver types</a>
     *              
     *              Setting own/utopia preconditioner, can be done as following: 
     *              \snippet tests/utopia_SolverTest.cpp PETScKSPSolver solve example1
     *              Detailed information about preconditioners, can be found in  \ref precondotioners.
     */
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend> 
    class KSPSolver {};


    template<typename Matrix, typename Vector>
    class KSPSolver<Matrix, Vector, PETSC> : virtual public IterativeSolver<Matrix, Vector>
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::Preconditioner<Vector> Preconditioner;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver; 

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "only works with petsc types");



        KSPSolver(  const Parameters params = Parameters(), 
                    const std::vector<std::string> ksp_types    = {"bcgs", "cg", "groppcg", "pipecg", "pipecgrr", "fcg", "pipefcg", "gmres", "pipefgmres",   "fgmres",   "lgmres",   "dgmres",   "pgmres", "tcqmr", "ibcgs",   "fbcgs",   "fbcgsr",   "bcgsl", "cgs", "tfqmr", "cr", "pipecr", "lsqr", "preonly", "qcg", "bicg", "minres", "symmlq", "lcd", "python", "gcr", "pipegcr", "tsirm", "cgls"},
                    const std::vector<std::string> pc_types     = {"jacobi","sor","lu","bjacobi","eisenstat","ilu","icc","asm","gasm","ksp","cholesky","pbjacobi","mat","hypre", "cp","bfbt","lsc","python","pfmg","syspfmg","redistribute","svd","gamg","bicgstabcusp","ainvcusp","bddc"}, 
                    const std::vector<std::string> pc_packages = {" "}):

                KSP_types(ksp_types),
                PC_types(pc_types),
                Solver_packages(pc_packages)

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
        IterativeSolver::set_parameters(params);

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


   /**
   * @brief      Sets the preconditioner.
   *
   * @param[in]  precond  The precondition
   */
    void set_preconditioner(const std::shared_ptr<Preconditioner> &precond)
   {
       precond_ = precond;
   }


    /**
    * @brief      Resets the preconditioner.
    */
    void reset_preconditioner()
    {
       precond_.reset(); 
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
    std::string pc_type() { return PC_type_;}; 

    /**
     * @brief      Returns ksp package type 
     */
    std::string ksp_type() { return KSP_type_; }; 


    /**
     * @brief      Returns type of solver package. 
     */
    std::string solver_package() { return solver_package_; }; 



public:

    /**
    * @brief        Solve function for all PETS KSP solvers. 
     *              It is also compatible with our own Utopia preconditioners.
    */
    bool apply(const Vector &b, Vector &x) override
    {
        PetscErrorCode ierr;

        Size localSize = local_size(b);

        if(empty(x) || local_size(x).get(0) != localSize.get(0)) 
        {
            x = local_zeros(localSize.get(0));
        }

        assert(b.size().get(0) == x.size().get(0));


        const Matrix A = *this->get_operator(); 


        IterativeSolver::init_solver("UTOPIA::PETSc KSP", {}); 
          
        KSPConvergedReason  reason;
        // PetscReal           r_norm;
        PetscInt            its; 
        MPI_Comm            comm; 
          

        ierr = PetscObjectGetComm((PetscObject)raw_type(A), &comm);
        ierr = KSPCreate(comm, &ksp);
        
        ierr = KSPSetOperators(ksp, raw_type(A), raw_type(A));

        set_ksp_options(ksp); 

        //  Set the user-defined routine for applying the preconditioner 
        if(precond_) 
        {
            PC pc; 
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCSHELL);  
            ierr = PCShellSetApply(pc, UtopiaPCApplyShell);
  
            auto shell_ptr = dynamic_cast<LinearSolver *>(precond_.get());
            ierr = PCShellSetContext(pc, shell_ptr);
            ierr = PCShellSetName(pc,"Utopia Preconditioner");
        }

        VecDuplicate(raw_type(x), &(ut_log.x_k_2));
        VecDuplicate(raw_type(x), &(ut_log.x_k_1));

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

    virtual void update(const std::shared_ptr<const Matrix> &op) override
    {
         IterativeSolver::update(op);

         if(precond_) {
            auto ls_ptr = dynamic_cast<LinearSolver *>(precond_.get());
            if(ls_ptr) {
                ls_ptr->update(op);    
            }
         }
     }

private: 

    virtual bool in_array(const std::string &value, const std::vector<std::string> &array)
    {
        return std::find(array.begin(), array.end(), value) != array.end();
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

        if(this->verbose())
            KSPMonitorSet(ksp, MyKSPMonitor, &ut_log, 0);
        
        ierr = KSPSetType(ksp, KSP_type_.c_str());

        if(!precond_) 
        {
            PC pc; 
            ierr = KSPGetPC(ksp, &pc);
            ierr = PCSetType(pc, PC_type_.c_str());
        }

        ierr = KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
        ierr = KSPSetTolerances(ksp, IterativeSolver::rtol(), IterativeSolver::atol(), PETSC_DEFAULT,  IterativeSolver::max_it());
    }


    std::shared_ptr<Preconditioner> precond_;

    std::string KSP_type_;                                  /*!< Choice of preconditioner types. */  
    const std::vector<std::string> KSP_types;              /*!< Valid options for direct solver types. */  

    std::string  PC_type_;                                  /*!< Choice of preconditioner types. */  
    const std::vector<std::string> PC_types;      

    std::string solver_package_;                            /*!< Choice of direct solver. */     
    const std::vector<std::string> Solver_packages;       /*!< Valid options for Solver packages types. */

protected:
    KSP                 ksp;
    UTOPIA_LOG          ut_log; 

    };
    
}



#endif //UTOPIA_PETSC_KSP_H
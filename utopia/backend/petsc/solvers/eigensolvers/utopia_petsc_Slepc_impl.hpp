#ifndef UTOPIA_PETSC_SLEPC_IMPL_HPP
#define UTOPIA_PETSC_SLEPC_IMPL_HPP

#include "utopia_petsc_Slepc.hpp"
#include <slepceps.h>

namespace utopia {

    template<typename Matrix, typename Vector>
    SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::SlepcSolver(
        const std::vector<std::string> problem_types,
        const std::vector<std::string> portions_of_spectrum,
        const std::vector<std::string> solver_types
        ) : 
        initialized_(false),
        solved_(false),
        problem_types_(problem_types),
        portions_of_spectrum_(portions_of_spectrum),
        solver_types_(solver_types)

    {
        problem_type_           = problem_types_.at(0);
        portion_of_spectrum_    = portions_of_spectrum_.at(0);
        solver_type_            = solver_types_.at(0);
    }

    template<typename Matrix, typename Vector>
    SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::~SlepcSolver()
    {
        if (initialized_)
            EPSDestroy(&eps_);
    }

    template<typename Matrix, typename Vector>
    SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL> * SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::clone() const 
    {
        return new SlepcSolver(*this);
    }

    template<typename Matrix, typename Vector>
    void SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::problem_type(const std::string & type)
    {
      problem_type_ = in_array(type, problem_types_) ? type : problem_types_.at(0);
    }

    template<typename Matrix, typename Vector>
    const std::string & SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::problem_type() const
    {
      return problem_type_;
    }

    template<typename Matrix, typename Vector>
    void SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::solver_type(const std::string & type)
    {
      solver_type_ = in_array(type, solver_types_) ? type : solver_types_.at(0);
    }

    template<typename Matrix, typename Vector>
    const std::string & SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::solver_type() const
    {
      return solver_type_;
    }

    template<typename Matrix, typename Vector>
    void SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::portion_of_spectrum(const std::string & type)
    {
      portion_of_spectrum_ = in_array(type, portions_of_spectrum_) ? type : portions_of_spectrum_.at(0);
    }

    template<typename Matrix, typename Vector>
    bool SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::solve(const Matrix & A)
    {
        MPI_Comm            comm;
        PetscObjectGetComm((PetscObject)raw_type(A), &comm);

        if (initialized_)
            reinitialize(comm);
        else
            initialize(comm);


        EPSSetOperators(eps_, raw_type(A), NULL);

        EPSSolve(eps_);

        EPSConvergedReason convergence_reason;
        EPSGetConvergedReason(eps_, &convergence_reason);

        // PetscInt its;
        // EPSGetIterationNumber(eps_, &its);
        // std::cout<<"it:  "<< its << "  \n";

        if(convergence_reason > 0)
            solved_ = true;
        else
            std::cerr << "[Error] EigenSolver did not converge." << std::endl;

        return true;
    }

    template<typename Matrix, typename Vector>
    bool SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::solve(const Matrix & A, const Matrix & B)
    {
        MPI_Comm            comm;
        PetscObjectGetComm((PetscObject)raw_type(A), &comm);

        if (initialized_)
            reinitialize(comm);
        else
            initialize(comm);


        EPSSetOperators(eps_, raw_type(A), raw_type(B));

        EPSSolve(eps_);

        EPSConvergedReason convergence_reason;
        EPSGetConvergedReason(eps_, &convergence_reason);

        // PetscInt its;
        // EPSGetIterationNumber(eps_, &its);
        // std::cout<<"it:  "<< its << "  \n";


        // std::cout<<"convergence_reason: " << convergence_reason << "  \n";

        if(convergence_reason > 0)
            solved_ = true;
        else
            std::cerr << "[Error] EigenSolver did not converge." << std::endl;

        return true;
    }

    template<typename Matrix, typename Vector>
    bool SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::print_eigenpairs()
    {

        if(!this->verbose())
            return false;

        if(!solved_) {
            std::cerr << "[Warning] You need to solve (or let the solver converge) before querying the eigenpairs.... \n";
            return false;
        }

        PetscInt nconv;

        PetscScalar         kr,ki;
        Vector              xr,xi;

        Mat A;

        EPSGetOperators(eps_, &A, NULL);

        MatCreateVecs(A, NULL, & raw_type(xr));
        MatCreateVecs(A, NULL, & raw_type(xi));


        EPSGetConverged(eps_, &nconv);
        std::string method = "EIGEN_PAIRS  \n Number of converged eigenpairs: " + std::to_string(nconv);
        PrintInfo::print_init(method, {"index", "real eigenvalue", "im eigenvalue"});

        if (nconv>0)
        {
          for (auto i=0; i<nconv; i++)
          {
            EPSGetEigenpair(eps_, i, &kr, &ki ,raw_type(xr), raw_type(xi));
            PrintInfo::print_iter_status(i, {kr, ki});
          }
          PetscPrintf(PETSC_COMM_WORLD,"\n");
        }

        return true;
    }


    template<typename Matrix, typename Vector>
    void SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::get_eigenpairs(const SizeType & i, Scalar & iegr, Scalar & eigi, Vector & vr, Vector & vi)
    {
        if(!solved_)
        {
            std::cerr << "[Warning] You need to solve (or let the solver converge) before querying the eigenpairs.... \n";
            return;
        }

        PetscInt nconv;

        Mat A;

        EPSGetOperators(eps_, &A, NULL);

        if(!empty(vr))
            VecDestroy(&raw_type(vr));
        if(!empty(vi))
            VecDestroy(&raw_type(vi));

        MatCreateVecs(A, NULL, & raw_type(vr));
        MatCreateVecs(A, NULL, & raw_type(vi));

        vr.implementation().set_initialized(true);
        vi.implementation().set_initialized(true);

        EPSGetConverged(eps_, &nconv);


        if (i < nconv) {
            EPSGetEigenpair(eps_, i, &iegr, &eigi ,raw_type(vr), raw_type(vi));
        }

    }

    template<typename Matrix, typename Vector>
    void SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::get_real_eigenpair(const SizeType & i, Scalar & iegr, Vector & vr)
    {
        if(i >= this->number_of_eigenvalues()) {
            std::cerr << "[Warning] Trying to access eigen-pair " << i << " which was not requested\n";
            return;
        }

        if(!solved_)
        {
            std::cerr << "[Warning] You need to solve (or let the solver converge) before querying the eigenpairs.... \n";
            // return;
        }

        PetscInt nconv;

        Mat A;
        EPSGetOperators(eps_, &A, nullptr);

        if(!empty(vr))
            VecDestroy(&raw_type(vr));


        MatCreateVecs(A, nullptr, &raw_type(vr));
        vr.implementation().set_initialized(true);

        EPSGetConverged(eps_, &nconv);


        // if (i < nconv)
        // 
        if (i >= nconv) {
            std::cerr << "[Warning] requested eigen pair was not resolved up to required accuracy\n";
            
        }

        EPSGetEigenpair(eps_, i, &iegr, nullptr, raw_type(vr), nullptr);
    }

    template<typename Matrix, typename Vector>
    void SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::initialize(const MPI_Comm & comm)
    {

        EPSCreate(comm, &eps_);

        if(problem_type_ == "hermitian")
            EPSSetProblemType(eps_, EPS_HEP);
        else if(problem_type_ == "generalized_hermitian")
            EPSSetProblemType(eps_, EPS_GHEP);
        else if(problem_type_ == "non_hermitian")
            EPSSetProblemType(eps_, EPS_NHEP);
        else if(problem_type_ == "generalized_non_hermitian")
            EPSSetProblemType(eps_, EPS_GNHEP);
        else if(problem_type_ == "generalized_hermitian_SPD_B")
            EPSSetProblemType(eps_, EPS_PGNHEP);
        else
            EPSSetProblemType(eps_, EPS_GHIEP);


        if(portion_of_spectrum_ == "largest_magnitude")
            EPSSetWhichEigenpairs(eps_, EPS_LARGEST_MAGNITUDE );
        else if(portion_of_spectrum_ == "smallest_magnitude")
            EPSSetWhichEigenpairs(eps_, EPS_SMALLEST_MAGNITUDE   );
        else if(portion_of_spectrum_ == "largest_real")
            EPSSetWhichEigenpairs(eps_,  EPS_LARGEST_REAL );
        else if(portion_of_spectrum_ == "smallest_real")
            EPSSetWhichEigenpairs(eps_, EPS_SMALLEST_REAL );
        else if(portion_of_spectrum_ == "largest_imaginary")
            EPSSetWhichEigenpairs(eps_, EPS_LARGEST_IMAGINARY );
        else if(portion_of_spectrum_ == "smallest_imaginary")
            EPSSetWhichEigenpairs(eps_, EPS_SMALLEST_IMAGINARY );
        else if(portion_of_spectrum_ == "closest_to_target")
            EPSSetWhichEigenpairs(eps_, EPS_TARGET_MAGNITUDE );
        else if(portion_of_spectrum_ == "closest_to_target_real")
            EPSSetWhichEigenpairs(eps_, EPS_TARGET_REAL );
        else if(portion_of_spectrum_ == "closest_to_target_imaginary")
            EPSSetWhichEigenpairs(eps_, EPS_TARGET_IMAGINARY );
        else if(portion_of_spectrum_ == "all_in_region")
            EPSSetWhichEigenpairs(eps_, EPS_ALL );


        EPSSetType(eps_, solver_type_.c_str());
        EPSSetTolerances(eps_, this->tol(), this->max_it());

        // could be done more sophisticated...
        EPSSetDimensions(eps_, this->number_of_eigenvalues(), PETSC_DEFAULT,  PETSC_DEFAULT);


        EPSSetFromOptions(eps_);

        initialized_ = true;
    }

    template<typename Matrix, typename Vector>
    void SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL>::reinitialize(const MPI_Comm & comm)
    {
        EPSDestroy(&eps_);
        initialize(comm);

        solved_ = false;
    }
}



#endif //UTOPIA_PETSC_SLEPC_IMPL_HPP
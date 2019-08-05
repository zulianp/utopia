#ifndef UTOPIA_SNES_SOLVER_HPP
#define UTOPIA_SNES_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
// #include "utopia_petsc.hpp"
#include "utopia_petsc_build_ksp.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"

#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscsnes.h>

namespace utopia
{
    // this should be potentially nonlinear solver ....
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class SNESSolver {};


    template<typename Matrix, typename Vector>
    class SNESSolver<Matrix, Vector, PETSC> : public NewtonBase<Matrix, Vector>,
                                              public NonLinearSmoother<Matrix, Vector>,
                                              public virtual VariableBoundSolverInterface<Vector>,
                                              public virtual Clonable //FIXME all non-linear solvers should be clonable
    {

    public:
        typedef UTOPIA_SCALAR(Vector)    Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;

        typedef typename NewtonBase<Matrix, Vector>::Solver  LinearSolver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>    Smoother;
        typedef utopia::NewtonBase<Matrix, Vector>           NonLinearSolver;
        typedef utopia::Function<Matrix, Vector>             Function;


        SNESSolver(const std::shared_ptr<LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
                   const std::vector<std::string> snes_types    = {"newtonls", "newtontr", "nrichardson", "ksponly", "vinewtonrsls", "vinewtonssls", "ngmres", "qn", "shell", "ngs", "ncg", "fas", "ms", "anderson"}):
        NonLinearSolver(linear_solver),
        SNES_types(snes_types),
        line_search_type_(SNESLINESEARCHBT) //SNESLINESEARCHBASIC)
        {
            SNES_type_ = SNES_types.at(0);
        }


        virtual SNESSolver * clone() const override
        {
            //FIXME make complete clone
            auto cloned_linear_solver = std::shared_ptr<LinearSolver>(this->linear_solver()->clone());
            auto cloned = utopia::make_unique<SNESSolver>(cloned_linear_solver);
            cloned->set_snes_type(SNES_type_);
            return cloned.release();
        }


        virtual ~SNESSolver()
        {

        }

        virtual void read(Input &in) override
        {
            NonLinearSolver::read(in);
            Smoother::read(in);

            std::string SNES_type_aux_;
            in.get("SNES_type", SNES_type_aux_);

            // checks if type is valid
            this->set_snes_type(SNES_type_aux_);
        }


        virtual void print_usage(std::ostream &os) const override
        {
            NonLinearSolver::print_usage(os);
            Smoother::print_usage(os);

            this->print_param_usage(os, "SNES_type", "string", "Type of Snes solver.", "newtonls");
        }


        virtual void set_snes_type(const std::string & type)
        {
            SNES_type_ = in_array(type, SNES_types) ? type : SNES_types.at(0);
        }


        virtual bool solve(Function &fun, Vector &x) override
        {
            using namespace utopia;

            SNES            snes;

            std::string method = "SNES_INTERFACE";
            this->init_solver(method, {"it", "|g|"});


            if (dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun) != nullptr)
            {
                PETSCUtopiaNonlinearFunction<Matrix, Vector> * fun_petsc = dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun);
                fun_petsc->getSNES(snes);
            }
            else
                setup_assembly_routines(snes, fun, x);

            set_snes_options(snes, this->atol(), this->rtol(), this->stol(), this->max_it());
            set_ksp(snes);
            set_variable_bounds(snes);


            SNESSolve(snes, NULL, raw_type(x));

            // exit solver
            PetscInt nonl_its;
            SNESGetIterationNumber(snes, &nonl_its);

            SNESConvergedReason reason;
            SNESGetConvergedReason(snes, &reason);

            this->exit_solver(nonl_its, reason);

            if (dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun) == nullptr)
                SNESDestroy(&snes);

            return true;
        }




        virtual bool smooth(Function & fun,  Vector &x, const Vector &rhs) override
        {
            using namespace utopia;

            SNES            snes;

            if (dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun) != nullptr)
            {
                PETSCUtopiaNonlinearFunction<Matrix, Vector> * fun_petsc = dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun);
                fun_petsc->getSNES(snes);
            }
            else
                setup_assembly_routines(snes, fun, x);


            set_snes_options(snes, 0.0, 0.0, 0.0, this->sweeps());
            set_ksp(snes);
            set_variable_bounds(snes);

            SNESSolve(snes, raw_type(rhs), raw_type(x));

            // needs to be reseted for use on other levels ...
            VecDestroy(&snes->vec_rhs);
            snes->vec_rhs =  NULL;

            if (dynamic_cast<PETSCUtopiaNonlinearFunction<Matrix, Vector> *>(&fun) == nullptr)
            {
                MatDestroy(&snes->jacobian);
                MatDestroy(&snes->jacobian_pre);

                SNESDestroy(&snes);
            }


            return true;

        }

        void set_line_search_type(const SNESLineSearchType ls_type)
        {
            line_search_type_ = ls_type;
        }

    protected:
        virtual void set_snes_options(SNES & snes,
                                      const Scalar & atol,
                                      const Scalar & rtol,
                                      const Scalar & stol,
                                      const SizeType & max_it)
        {
            PetscErrorCode ierr; UTOPIA_UNUSED(ierr);

            SNESMonitorCancel(snes);
            SNESSetFromOptions(snes);


            if(this->verbose()) {
                SNESMonitorSet(
                               snes,
                               [](SNES snes, PetscInt iter, PetscReal res, void*) -> PetscErrorCode
                               {
                                   if(mpi_world_rank() == 0)
                                        std::cout<<iter << "       "<< res << "      \n";

                                   return 0;
                               },
                               nullptr,
                               nullptr);
            }

            ierr = SNESSetType(snes, SNES_type_.c_str());
            ierr = SNESSetTolerances(snes, atol, rtol, stol, max_it, PETSC_DEFAULT);


            SNESLineSearch linesearch;
            SNESGetLineSearch(snes, &linesearch);

            SNESLineSearchSetFromOptions(linesearch);
            SNESLineSearchSetType(linesearch, line_search_type_);

#if !UTOPIA_PETSC_VERSION_LESS_THAN(3,8,0)

            PetscReal damping = 1.;
            SNESLineSearchGetDamping(linesearch, &damping);
            SNESLineSearchSetTolerances(linesearch, PETSC_DEFAULT, PETSC_DEFAULT, this->rtol(), this->atol(), PETSC_DEFAULT, 20);

            // if(std::string(SNESLINESEARCHBASIC) == line_search_type_ && std::abs(damping - 1.) < 1e-16) {
            //       SNESLineSearchSetComputeNorms(linesearch, PETSC_FALSE);
            // }
#endif
        }



        virtual void set_ksp(SNES & snes)
        {

            if(!snes->usesksp)
                return;

            KSP            ksp;
            SNESGetKSP(snes,&ksp);

            if (dynamic_cast<KSPSolver<Matrix, Vector>*>(this->linear_solver_.get()) != nullptr)
            {
                auto utopia_ksp = dynamic_cast<KSPSolver<Matrix, Vector> *>(this->linear_solver_.get());
                utopia_ksp->set_ksp_options(ksp);
                utopia_ksp->attach_preconditioner(ksp);
            }
            else
            {
                if(!this->linear_solver_)
                {
                    std::cout<<"utopia::SNES:: linear solver missing, setting to DEFAULT... \n";
                    const auto utopia_ksp = std::make_shared<KSPSolver<Matrix, Vector> >();
                    this->set_linear_solver(utopia_ksp);
                }

                build_ksp(this->linear_solver_, ksp);
            }
        }


    private:

        void set_variable_bounds(SNES &snes)
        {
            if(this->has_bound()) {
                this->fill_empty_bounds();

                SNESVISetVariableBounds(
                    snes,
                    raw_type(this->get_lower_bound()),
                    raw_type(this->get_upper_bound())
                );
            }
        }

        // this function is expensive - wrt to convert
        // however, in utopia is not ready for pure wrap implementation
        void setup_assembly_routines(SNES & snes, Function & fun, const Vector &x)
        {
            MPI_Comm        comm;
            PetscObjectGetComm((PetscObject)raw_type(x), &comm);

            // residual
            Vector residual = local_zeros(local_size(x));

            SNESCreate(comm, &snes);

            // energy
            SNESSetObjective( snes,
                             // FormObjective,
                             [](SNES /*snes*/, Vec x, PetscReal * energy, void * ctx) -> PetscErrorCode
                             {
                                 Function * fun = static_cast<Function*>(ctx);
                                 Vector x_ut;

                                 utopia::convert(x, x_ut);
                                 fun->value(x_ut, *energy);

                                 return 0;
                             },

                             &fun);

            // gradient
            SNESSetFunction( snes,
                            raw_type(residual),
                            // FormGradient,
                            [](SNES snes, Vec x, Vec res, void *ctx)-> PetscErrorCode
                            {
                                Function * fun = static_cast<Function *>(ctx);

                                Vector x_ut, res_ut;
                                utopia::convert(x, x_ut);
                                fun->gradient(x_ut, res_ut);
                                utopia::convert(res_ut, res);

                                return 0;
                            },
                            &fun);

            // hessian
            SNESSetJacobian( snes,  snes->jacobian,  snes->jacobian_pre,
                            // FormHessian,
                            [](SNES snes, Vec x, Mat jac, Mat prec, void *ctx)-> PetscErrorCode
                            {
                                Function * fun = static_cast<Function *>(ctx);

                                Vector x_ut;
                                utopia::convert(x, x_ut);

                                Matrix jac_ut, jac_ut_prec;

                                bool flg = fun->hessian(x_ut, jac_ut, jac_ut_prec);

                                if(!flg)
                                {
                                    fun->hessian(x_ut, jac_ut);
                                    convert(jac_ut, snes->jacobian_pre);
                                }
                                else
                                    convert(jac_ut_prec, snes->jacobian_pre);

                                convert(jac_ut, snes->jacobian);

                                return 0;
                            },
                            &fun);
        }


    protected:
        std::string SNES_type_;                                  /*!< Choice of snes types. */
        const std::vector<std::string> SNES_types;              /*!< Valid options for SNES solver types. */
        SNESLineSearchType line_search_type_;


    };
}

#endif // UTOPIA_SNES_SOLVER_HPP

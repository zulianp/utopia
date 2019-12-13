#ifndef UTOPIA_BRATU2D_HPP
#define UTOPIA_BRATU2D_HPP

#include "utopia.hpp"
#include "utopia_TestFunctions.hpp"



namespace utopia
{
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Bratu2D  { }; 
}

#ifdef  WITH_PETSC
#warning "this code should go in the petsc backend folder"

#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscmatlab.h>
#include <petsc/private/snesimpl.h> /* For SNES_Solve event */

namespace utopia
{
    typedef struct ParamsBratu2D ParamsBratu2D;
    struct ParamsBratu2D 
    {
      PetscReal lambda;          /* test problem parameter */
      PetscErrorCode (*mms_solution)(ParamsBratu2D*,const DMDACoor2d*,PetscScalar*);
    };

    PetscErrorCode Bratu2DFormObjectiveLocal(DMDALocalInfo*,PetscScalar**,PetscReal*,ParamsBratu2D*);
    PetscErrorCode Bratu2DFormFunctionLocal(DMDALocalInfo*,PetscScalar**,PetscScalar**,ParamsBratu2D*);
    PetscErrorCode Bratu2DFormJacobianLocal(DMDALocalInfo*,PetscScalar**,Mat,Mat,ParamsBratu2D*);
    
    PetscErrorCode Bratu2DFormInitialGuess(DM,ParamsBratu2D*,Vec);
    PetscErrorCode Bratu2DMMSSolution(ParamsBratu2D*,const DMDACoor2d*,PetscScalar*);

    PetscErrorCode Bratu2DFormBCData(DM da,ParamsBratu2D *user,Vec BC_marker, Vec BC_flag);


    template<typename Matrix, typename Vector>
    class Bratu2D<Matrix, Vector, PETSC> final: virtual public UnconstrainedExtendedTestFunction<Matrix, Vector>, virtual public ConstrainedExtendedTestFunction<Matrix, Vector>
    {
        public:
            typedef UTOPIA_SIZE_TYPE(PetscVector) SizeType;
            typedef UTOPIA_SCALAR(PetscVector) Scalar;


        Bratu2D(const SizeType & n,
                const Scalar & lambda = 5.0):
                n_(n), 
                setup_(false)
        {
            application_context_.lambda  = (lambda >= 0 && lambda <= 6.8) ? lambda : 3.4; 
            application_context_.mms_solution   = Bratu2DMMSSolution; 

            this->create_DM();
            this->setup_SNES();
            this->setup_application_context(); 
            setup_ = true;

            PetscInt n_local_; 
            VecGetLocalSize(snes_->vec_sol, &n_local_);
            this->constraints_ = make_box_constaints(std::make_shared<Vector>(local_values(n_local_, -9e9)),
                                                     std::make_shared<Vector>(local_values(n_local_, 0.45)));    


            // TD::fix this
            exact_sol_  = zeros(1, 0); 
        }     

        Bratu2D(const DM  & dm, const Scalar & lambda = 5.0):
                setup_(false)
        {
            application_context_.lambda  = (lambda >= 0 && lambda <= 6.8) ? lambda : 3.4; 
            application_context_.mms_solution   = Bratu2DMMSSolution; 
            // application_context_.mms_forcing    = Bratu2DMMSForcing;

            da_ = dm; 
            // necessary to provide reasonable global dimension 
            DMDAGetInfo(da_, 0, &n_, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

            this->setup_SNES();
            this->setup_application_context(); 
            setup_ = true;

            // TODO::find out exact solution - should be possible to compute
            exact_sol_  = zeros(1, 0); 
        }     

        ~Bratu2D()
        {
            if(setup_)
            {
                DMDestroy(&da_);
                SNESDestroy(&snes_);
            }
        }

        virtual bool gradient(const Vector &x, Vector &g) const override
        {
            // initialization of gradient vector...
            if(empty(g)){
                g  = local_zeros(local_size(x));;
            }

            SNESComputeFunction(snes_, raw_type(x), raw_type(g));

            // remove BC from gradient... 

            return true;
        }
            
        virtual bool hessian(const Vector &x, Matrix &hessian) const override
        {
            
            SNESComputeJacobian(snes_, raw_type(x), snes_->jacobian, snes_->jacobian);

            // yes, wrap would be nicer, but lets not use it ...
            convert(snes_->jacobian, hessian); 

            // disp(hessian); 
            // write("hessian.m", hessian); 
            
            return true;
        }

        virtual bool value(const Vector &x, typename Vector::Scalar &result) const override
        {
            SNESComputeObjective(snes_, raw_type(x), &result);
            return true;
        }


        void output_to_VTK(const Vector & x, const std::string file_name = "Bratu2D.vtk")
        {
          PetscViewer       viewer;

          PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
          PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_VTK);
          DMView(da_, viewer);
          PetscObjectSetName((PetscObject)raw_type(x), "x");
          VecView(raw_type(x), viewer);
          PetscViewerDestroy(&viewer);
        }


        virtual Vector initial_guess() const override
        {   
            Vector x_utopia; 
            convert(snes_->vec_sol, x_utopia); 
            return x_utopia; 
        }
        
        virtual const Vector & exact_sol() const override
        {
            return exact_sol_; 
        }
        

        virtual Scalar min_function_value() const override
        {   
            // depends on the solution to which we converged to 
            return -1.012; 
        }

        virtual std::string name() const override
        {
            return "Bratu2D";
        }
        
        virtual SizeType dim() const override
        {
            return n_*n_; 
        }

        virtual bool exact_sol_known() const override
        {
            return true;
        }

        virtual bool parallel() const override
        {
            return true;
        }


        virtual void lambda(const Scalar & lambda)
        {
            application_context_.lambda = lambda;
        }

        virtual Scalar lambda() const
        {
            return application_context_.lambda; 
        }


    private:
        void create_DM()
        {
            DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,DMDA_STENCIL_STAR, n_, n_,PETSC_DECIDE,PETSC_DECIDE,1,1,NULL,NULL, &da_);
            DMSetUp(da_);
            DMDASetUniformCoordinates(da_, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
        }

        void setup_SNES()
        {
            SNESCreate(PETSC_COMM_WORLD, &snes_); 
            SNESSetFromOptions(snes_);
            SNESSetDM(snes_, da_);
            DMDASNESSetObjectiveLocal(da_,(DMDASNESObjective)Bratu2DFormObjectiveLocal,&application_context_);            
            DMDASNESSetFunctionLocal(da_,INSERT_VALUES,(DMDASNESFunction)Bratu2DFormFunctionLocal,&application_context_);
            DMDASNESSetJacobianLocal(da_,(DMDASNESJacobian)Bratu2DFormJacobianLocal,&application_context_);
            
            // preallocate matrices/vectors 
            DMCreateMatrix(da_, &snes_->jacobian);
            DMCreateMatrix(da_, &snes_->jacobian_pre);
            DMCreateGlobalVector(da_, &snes_->vec_sol);
            Bratu2DFormInitialGuess(da_, &application_context_, snes_->vec_sol);
        }

        void setup_application_context()
        {
            DMSetApplicationContext(da_, &application_context_);

            PetscInt n_loc; 
            VecGetLocalSize(snes_->vec_sol, &n_loc); 
            Vector bc_markers = local_values(n_loc, 0.0);
            Vector bc_values  = local_values(n_loc, 0.0); 

            Bratu2DFormBCData(da_, &application_context_, raw_type(bc_markers), raw_type(bc_values)); 
            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, bc_values);
        }


    private:
        SizeType n_;  // global size
        bool setup_; 

        ParamsBratu2D application_context_; 
        DM da_;
        SNES snes_; 

        Vector exact_sol_; 

    };
}

#endif //WITH_PETSC
#endif

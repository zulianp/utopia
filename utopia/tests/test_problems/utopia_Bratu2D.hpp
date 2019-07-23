#include "utopia.hpp"
#include "utopia_UnconstrainedTestFunction.hpp"

namespace utopia
{
    template<class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class Bratu2D  { }; 
}


#ifdef  WITH_PETSC
#include <petscdm.h>
#include <petscdmda.h>
#include <petscsnes.h>
#include <petscmatlab.h>
#include <petsc/private/snesimpl.h> /* For SNES_Solve event */

namespace utopia
{
    typedef struct AppCtxBratu2D AppCtxBratu2D;
    struct AppCtxBratu2D 
    {
      PetscReal lambda;          /* test problem parameter */
      PetscErrorCode (*mms_solution)(AppCtxBratu2D*,const DMDACoor2d*,PetscScalar*);
      PetscErrorCode (*mms_forcing)(AppCtxBratu2D*,const DMDACoor2d*,PetscScalar*);
    };

    PetscErrorCode FormObjectiveLocal(DMDALocalInfo*,PetscScalar**,PetscReal*,AppCtxBratu2D*);
    PetscErrorCode FormFunctionLocal(DMDALocalInfo*,PetscScalar**,PetscScalar**,AppCtxBratu2D*);
    PetscErrorCode FormJacobianLocal(DMDALocalInfo*,PetscScalar**,Mat,Mat,AppCtxBratu2D*);
    
    PetscErrorCode FormInitialGuess(DM,AppCtxBratu2D*,Vec);
    PetscErrorCode FormExactSolution(DM,AppCtxBratu2D*,Vec);
    PetscErrorCode MMSSolution(AppCtxBratu2D*,const DMDACoor2d*,PetscScalar*);
    PetscErrorCode MMSForcing(AppCtxBratu2D*,const DMDACoor2d*,PetscScalar*);    

    PetscErrorCode FormBCData(DM da,AppCtxBratu2D *user,Vec BC_marker, Vec BC_flag);


    template<typename Matrix, typename Vector>
    class Bratu2D<Matrix, Vector, PETSC> : public UnconstrainedExtendedTestFunction<Matrix, Vector>
    {
        public:
            typedef UTOPIA_SIZE_TYPE(DVectord) SizeType;
            typedef UTOPIA_SCALAR(DVectord) Scalar;


        Bratu2D(const SizeType & n,
                const Scalar & lambda = 0.1):
                n_(n), 
                setup_(false)
        {
            application_context_.lambda  = (lambda > 0 && lambda < 6.8) ? lambda : 3.4; 
            application_context_.mms_solution = MMSSolution; 
            application_context_.mms_forcing = MMSForcing;

            this->create_DM();
            this->setup_SNES();
            this->setup_application_context(); 
            setup_ = true;
        }     

        Bratu2D(const DM  & dm,
                const AppCtxBratu2D & context):
                setup_(false)
        {
            application_context_ = context; 
            da_ = dm; 
            DMDAGetInfo(da_, 0, &n_, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

            this->setup_SNES();
            this->setup_application_context(); 
            setup_ = true;
        }     

        ~Bratu2D()
        {
            if(setup_)
            {
                DMDestroy(&da_);
                SNESDestroy(&snes_);
            }
        }

            virtual bool gradient_no_rhs(const Vector &x, Vector &g) const override
            {
                // initialization of gradient vector...
                if(empty(g)){
                    g  = local_zeros(local_size(x));;
                }

                SNESComputeFunction(snes_, raw_type(x), raw_type(g));

                return true;
            }
                
            virtual bool hessian(const Vector &x, Matrix &hessian) const override
            {
                SNESComputeJacobian(snes_, raw_type(x), snes_->jacobian,  snes_->jacobian);
                wrap(snes_->jacobian, hessian);
                return true;
            }

            virtual bool value(const Vector &x, typename Vector::Scalar &result) const override
            {
                SNESComputeObjective(snes_, raw_type(x), &result);
                return true;
            }


        void output_to_VTK(const Vector & x, const std::string file_name = "Bratu2D.vtk")
        {
          PetscErrorCode ierr;
          PetscViewer       viewer;

          PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
          PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_VTK);
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
            Vector empty; 
            return empty; 
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
            DMDASNESSetFunctionLocal(da_,INSERT_VALUES,(DMDASNESFunction)FormFunctionLocal,&application_context_);
            DMDASNESSetJacobianLocal(da_,(DMDASNESJacobian)FormJacobianLocal,&application_context_);
            DMDASNESSetObjectiveLocal(da_,(DMDASNESObjective)FormObjectiveLocal,&application_context_);
            
            // preallocate vectors 
            DMCreateMatrix(da_, &snes_->jacobian);
            DMCreateMatrix(da_, &snes_->jacobian_pre);
            DMCreateGlobalVector(da_, &snes_->vec_sol);
            FormInitialGuess(da_, &application_context_, snes_->vec_sol);
        }

        void setup_application_context()
        {
            DMSetApplicationContext(da_, &application_context_);

            PetscInt n_loc; 
            VecGetLocalSize(snes_->vec_sol, &n_loc); 
            Vector bc_markers = local_values(n_loc, 0.0);
            Vector bc_values  = local_values(n_loc, 0.0); 

            FormBCData(da_, &application_context_, raw_type(bc_markers), raw_type(bc_values)); 
            ExtendedFunction<Matrix, Vector>::set_equality_constrains(bc_markers, bc_values);
        }

    private:
        SizeType n_;  // global size
        bool setup_; 

        AppCtxBratu2D application_context_; 
        DM da_;
        SNES snes_; 

    };
}

#endif //WITH_PETSC

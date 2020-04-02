#ifndef PETSC_BASED_UTOPIA_NONLINEAR_FUNCTION_HPP
#define PETSC_BASED_UTOPIA_NONLINEAR_FUNCTION_HPP

#include "utopia_Core.hpp"
#include "utopia_ExtendedFunction.hpp"

#include "utopia_petsc_Types.hpp"

#include <petscsnes.h>
#include <petsc/private/snesimpl.h>

namespace utopia
{
    template<class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class PETSCUtopiaNonlinearFunction {};


    template<class Matrix, class Vector>
    class PETSCUtopiaNonlinearFunction<Matrix, Vector, PETSC> : public ExtendedFunction<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)    Scalar;

        public:
            // PETSCUtopiaNonlinearFunction(SNES snes, const Vector & x_init = local_zeros(1), const Vector & bc_marker = local_zeros(1), const Vector & rhs = local_zeros(1)) :
            PETSCUtopiaNonlinearFunction(SNES snes, const Vector & x_init = Vector(), const Vector & bc_marker = Vector(), const Vector & rhs = Vector()) :
                ExtendedFunction<Matrix, Vector>(x_init, bc_marker, rhs), snes_(snes)
            {}

            virtual bool gradient(const Vector &x, Vector &g) const override
            {
                // initialization of gradient vector...
                if(empty(g)){
                    g.zeros(layout(x));
                }

                SNESComputeFunction(snes_, raw_type(x), raw_type(g));

                return true;
            }

            virtual bool hessian(const Vector &x, Matrix &hessian) const override
            {
                SNESComputeJacobian(snes_, raw_type(x), snes_->jacobian,  snes_->jacobian_pre);
                wrap(snes_->jacobian, hessian);
                return true;
            }

            virtual bool value(const Vector &x, typename Vector::Scalar &result) const override
            {
                // hack to have fresh energy (MOOSE post-processor does things in strange way )
                Vector grad = 0 * x;
                this->gradient(x, grad);


                DM dm;
                DMSNES         sdm;

                SNESGetDM(snes_,&dm);
                DMGetDMSNES(dm,&sdm);
                if (sdm->ops->computeobjective)
                    SNESComputeObjective(snes_, raw_type(x), &result);
                else
                    result = 0.5 * norm2(grad) * norm2(grad);


                return true;
            }

            virtual void  getSNES(SNES &snes)
            {
                snes = snes_;
            }


        private:

            SNES snes_;
        };

    }


#endif  //PETSC_BASED_UTOPIA_NONLINEAR_FUNCTION_HPP

/*
* @Author: alenakopanicakova
* @Date:   2017-04-17
* @Last Modified by:   Alena Kopanicakova
* @Last Modified time: 2017-07-03
*/

#ifndef UTOPIA_PETSC_NONLINEAR_SMOOTHER_HPP
#define UTOPIA_PETSC_NONLINEAR_SMOOTHER_HPP

#include "utopia_Smoother.hpp"
#include "utopia_Core.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_NonLinearSmoother.hpp"

#include "utopia_petsc_SNESFunction.hpp"
#include "utopia_petsc_SNES.hpp"

#include <petsc/private/snesimpl.h>
#include "petscsnes.h"  

#include "utopia_NonlinearSolverInterfaces.hpp"


// TODO:: check nonlinear preconditioners ... 
namespace utopia 
{

    template<class Matrix, class Vector>
    class NonLinearGaussSeidel<Matrix, Vector, PETSC> : public SNESSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

        typedef utopia::SNESSolver<Matrix, Vector>                  SNESSolver;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver    LinearSolver;


        public:
        NonLinearGaussSeidel(   const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(), 
                                const Parameters params = Parameters()) 
                                : SNESSolver(linear_solver, params)
        { 
            set_parameters(params); 
            this->set_snes_type("ngs"); 

        }

        virtual void set_parameters(const Parameters params) override
        {
            SNESSolver::set_parameters(params); 
        }

    protected: 
        virtual void set_snes_options(SNES & snes,  const Scalar & atol     = SNESSolver::atol(), 
                                                    const Scalar & rtol     = SNESSolver::rtol(), 
                                                    const Scalar & stol     = SNESSolver::stol(), 
                                                    const SizeType & max_it = SNESSolver::max_it()) override 
        {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it); 

            // we need to allocate hessian for coloring computation

            PetscBool assembled; 
            MatAssembled(snes->jacobian_pre, &assembled); 

            if(!assembled)
                SNESComputeJacobian(snes, snes->vec_sol, snes->jacobian,  snes->jacobian_pre);

            SNESLineSearch linesearch; 
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch, SNESLINESEARCHBASIC); 
        }

    };



    template<class Matrix, class Vector>
    class NonLinearConjugateGradient<Matrix, Vector, PETSC> : public SNESSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)                           Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)                        SizeType;

        typedef utopia::SNESSolver<Matrix, Vector>                  SNESSolver;
        typedef typename NonLinearSolver<Matrix, Vector>::Solver    LinearSolver;


        public:
        NonLinearConjugateGradient( const std::shared_ptr <LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(), 
                                    const Parameters params = Parameters(), 
                                    const std::vector<std::string> update_types    = {"FR", "PRP", "HS", "DY", "CD"}) 
                                    :   SNESSolver(linear_solver, params), 
                                        update_types(update_types)
        { 
            set_parameters(params); 
            this->set_snes_type("ncg"); 
            update_type_ = update_types.at(0); 
        }

        virtual void set_parameters(const Parameters params) override
        {
            SNESSolver::set_parameters(params); 
        }


        void set_update_type(const std::string & update_type )
        {
            update_type_ = in_array(update_type, update_types) ? update_type : update_types.at(0);
        }

        std::string & get_update_type()
        {
            return update_type_; 
        }


    protected: 
        virtual void set_snes_options(SNES & snes,  const Scalar & atol     = SNESSolver::atol(), 
                                                    const Scalar & rtol     = SNESSolver::rtol(), 
                                                    const Scalar & stol     = SNESSolver::stol(), 
                                                    const SizeType & max_it = SNESSolver::max_it()) override 
        {
            SNESSolver::set_snes_options(snes, atol, rtol, stol, max_it); 

            if(update_type_ == "PRP")
                SNESNCGSetType(snes, SNES_NCG_PRP); 
            else if(update_type_ == "HS")
                SNESNCGSetType(snes, SNES_NCG_HS ); 
            else if(update_type_ == "DY")
                SNESNCGSetType(snes, SNES_NCG_DY); 
            else if(update_type_ == "CD")
                SNESNCGSetType(snes, SNES_NCG_CD ); 
            else
                SNESNCGSetType(snes, SNES_NCG_FR ); 

            // error oriented LS seems to work the best ... 
            SNESLineSearch linesearch; 
            SNESGetLineSearch(snes, &linesearch);
            SNESLineSearchSetType(linesearch,   SNESLINESEARCHCP   ); 
        }


    private: 
        const std::vector<std::string> update_types; 
        std::string update_type_; 

    };












}

#endif //UTOPIA_PETSC_NONLINEAR_SMOOTHER_HPP


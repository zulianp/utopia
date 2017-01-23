/*
* @Author: alenakopanicakova
* @Date:   2016-05-09
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-05-27
*/

#ifndef FENICS_UTOPIA_HPP
#define FENICS_UTOPIA_HPP

#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include "Form.h"
#include <utopia.hpp>

namespace utopia 
{

template<class Matrix, class Vector> 
class FenicsUtopiaFunction : public utopia::Function<Matrix, Vector> 
{



    public:
        FenicsUtopiaFunction(   std::shared_ptr<dolfin::Function>   u,
                                std::shared_ptr<const dolfin::Form> E, 
                                std::shared_ptr<const dolfin::Form> g, 
                                std::shared_ptr<const dolfin::Form> H,
                                std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs):
            _u(u),
            _E(E),  
            _g(g),
            _H(H), 
            _bcs(bcs)
        {

        }

        bool gradient(const Vector &x, Vector &g) const override
        {

            // this is changing u in lin and bil. forms 
            dolfin::PETScVector x_wrap(utopia::raw_type(x)); 
            x_wrap.update_ghost_values();
            (*_u->vector()) = x_wrap;

            dolfin::PETScVector b;
            dolfin::SystemAssembler assembler(_H, _g, _bcs);
            assembler.assemble(b, *_u->vector());

            Vec up = b.vec(); 
            convert(up, g); 

            return true; 
        }



        bool hessian(const Vector &x, Matrix &H) const override
        {

            // this is changing u in lin and bil. forms 
            dolfin::PETScVector x_wrap(utopia::raw_type(x)); 
            x_wrap.update_ghost_values();

            dolfin::PETScMatrix A;
            dolfin::SystemAssembler assembler(_H, _g, _bcs);
            assembler.assemble(A);

            Mat Ap = A.mat(); 
            // convert(Ap, H); 
            H = sparse_mref(Ap);
            
            return true; 
        }



        // energy 
        bool value(const Vector &x, typename Vector::Scalar &f) const override 
        {
            // this is changing u in lin and bil. forms 
            dolfin::PETScVector x_wrap(utopia::raw_type(x)); 
            (*_u->vector()) = x_wrap;

            dolfin::Scalar global_energy;
            dolfin::assemble(global_energy, *_E);             
            f  = global_energy.get_scalar_value(); 
            
            return true; 
        }


    private:    
        std::shared_ptr<dolfin::Function>   _u;                             // solution
        std::shared_ptr<const dolfin::Form> _E;                             // energy 
        std::shared_ptr<const dolfin::Form> _g;                             // gradient
        std::shared_ptr<const dolfin::Form> _H;                             // hessian
        std::vector<std::shared_ptr<const dolfin::DirichletBC>> _bcs;       // BC 

    };

}




#endif  //FENICS_UTOPIA_HPP
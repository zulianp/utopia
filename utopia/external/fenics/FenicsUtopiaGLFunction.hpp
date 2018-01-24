/*
* @Author: alenakopanicakova
* @Date:   2016-05-09
* @Last Modified by:   alenakopanicakova
* @Last Modified time: 2016-06-06
*/

// .................................................................... WORK IN PROGRESS ...................................................................

#ifndef FENICS_UTOPIA_GLOBAL_LOCAL_HPP
#define FENICS_UTOPIA_GLOBAL_LOCAL_HPP

#include <dolfin/function/Function.h>
#include <dolfin/la/GenericVector.h>
#include "Form.h"
#include <utopia.hpp>
#include "FenicsUtopiaDecomposition.hpp"

namespace utopia 
{

template <class GlobalMatrix, class GlobalVector, class LocalMatrix, class LocalVector>
class FenicsUtopiaGLFunction : public utopia::GLFunction<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> 
{

    typedef utopia::FenicsDecomposition<GlobalMatrix, GlobalVector, LocalMatrix, LocalVector> FenicsDecomposition;

    public:
        FenicsUtopiaGLFunction( std::shared_ptr<dolfin::Function>   u,
                                std::shared_ptr<const dolfin::Form> E, 
                                std::shared_ptr<const dolfin::Form> g, 
                                std::shared_ptr<const dolfin::Form> H,
                                std::vector<std::shared_ptr<const dolfin::DirichletBC>> bcs,
                                const std::shared_ptr<FenicsDecomposition> &decomposition = std::shared_ptr<FenicsDecomposition>()):
            _u(u),
            _E(E),  
            _g(g),
            _H(H), 
            _bcs(bcs), 
            _decomposition(decomposition)

        {

        }


        bool init(const GlobalVector &x) const override
        {
            VecGetLocalSize(utopia::raw_type(x), &local_dim);
            this->_decomposition->init(x); 
            return true; 
        }

        bool gradient(const GlobalVector &x, GlobalVector &g) const override
        {

            dolfin::PetscVector x_wrap(utopia::raw_type(x)); 
            x_wrap.update_ghost_values();
            (*_u->vector()) = x_wrap;

            dolfin::PetscVector b;
            dolfin::SystemAssembler assembler(_H, _g, _bcs);
            assembler.assemble(b, *_u->vector());

            Vec up = b.vec(); 
            convert(up, g); 
            return true; 
        }



        bool hessian(const GlobalVector &x, GlobalMatrix &H) const override
        {

            // this is changing u in lin and bil. forms 
            dolfin::PetscVector x_wrap(utopia::raw_type(x)); 
            (*_u->vector()) = x_wrap;

            dolfin::PetscMatrix A;
            dolfin::SystemAssembler assembler(_H, _g, _bcs);
            assembler.assemble(A);

            Mat Ap = A.mat(); 
            convert(Ap, H); 
            return true; 
        }



        // energy 
        bool value(const GlobalVector &x, typename GlobalVector::Scalar &f) const override 
        {
            // this is changing u in lin and bil. forms 
            dolfin::PetscVector x_wrap(utopia::raw_type(x)); 
            (*_u->vector()) = x_wrap;

            dolfin::Scalar global_energy;
            dolfin::assemble(global_energy, *_E);             
            f  = global_energy.get_scalar_value(); 
            
            return true; 
        }


// -----------------------------------------------------------------------------------------------------------------------
// --------------------------------------------- local evaluations  ------------------------------------------------------
//                                  subscript k stands for quantity on k-th subset     
// -----------------------------------------------------------------------------------------------------------------------

        // this is so far working with global u update .. 
        //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CHANGE THIS  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        bool local_value(const LocalVector & x_k, typename LocalVector::Scalar & f_k) const  override 
        {   
            GlobalVector x; 
            this->interpolate(x_k, x);

            // this is changing u in lin and bil. forms 
            dolfin::PetscVector x_wrap(utopia::raw_type(x)); 
            (*_u->vector()) = x_wrap;


            dolfin::TensorLayout tl(MPI_COMM_SELF,{}, 0, {}, {});
            dolfin::Scalar local_value;
            local_value.init(tl);
            dolfin::assemble(local_value, *_E);

            // checking for NAN's
            if(local_value.get_scalar_value() != local_value.get_scalar_value())
            {

            }
            else
            {
                f_k = local_value.get_scalar_value(); 
            }

            return true;
        }


        bool local_gradient(const LocalVector & x_k, LocalVector & g_k ) const override 
        {

            GlobalVector g, x; 

            this->interpolate(g_k, g);
            this->interpolate(x_k, x);

            this->gradient(x, g); 
            this->restrict(g, g_k); 

            // // utopia::disp(g); 
            // // exit(0); 
            // if(mpi_rank == 1 )
            // {
            //     std::cout<<"local one ---------------------- \n"; 
            //     utopia::disp(g_k);             
            // }
            // exit(0); 
            
            // Vec l; 
            // VecGhostGetLocalForm(utopia::raw_type(g), &l); 

            // GlobalVector lu; 
            // convert(l, lu); 
            // //disp(g); 
            // std::cout<<"size lu: "<< lu.size().get(0)<< "  \n"; 
            // exit(0);

            return true;
        }



        bool local_hessian(const LocalVector &x_k, LocalMatrix &H_k) const override 
        {

            // TODO: use local block diag ... 

            GlobalVector x; 
            GlobalMatrix J, J_u;


            this->interpolate(H_k, J);
            this->interpolate(x_k, x);
            
            this->hessian(x, J); 
            Mat &J_k = utopia::raw_type(J);

            Mat _a; 

            // to get diag blocks - for local matrices 
            MatGetDiagonalBlock(J_k, &_a); 
            utopia::convert(_a, J_u);


            Range rr = rowRange(J_u);
            const SizeType rb = rr.begin();
            const SizeType local_extent_r = rr.extent();

            // assign from local Petsc matrixes into BLAS/utopia mat 
            {
                Read<GlobalMatrix> w(J_u);

                for (SizeType i = 0; i < local_extent_r  ; ++i)
                {
                    for (SizeType j = 0; j < local_extent_r  ; ++j)
                    {
                        H_k.set(i , j , J_u.get(rb + i, rb + j));
                    }
                }
            }


            // // utopia::disp(J); 
            // // exit(0); 
            // if(mpi_rank == 0)
            // {
            //     std::cout<<"local one ---------------------- \n"; 
            //     utopia::disp(H_k);             
            // }
            //exit(0); 


            return true;

        }




//-----------------------------------------  DECOMPOSITION ------------------------------------------

        // projecting local vector into global 
        bool interpolate(const LocalVector & x_k, GlobalVector & x) const override 
        {
            this->_decomposition->interpolate(x_k, x);
            return true;
        }

        // interpolating local sub-matrixes into global matrix 
        bool interpolate(const LocalMatrix & M_k, GlobalMatrix & M) const override
        {
            this->_decomposition->interpolate(M_k, M);
            return true; 
        }


        // projecting global vector into local 
        bool restrict(const GlobalVector &x, LocalVector &x_k) const override 
        {
            this->_decomposition->restrict(x, x_k);
            return true;
        }

        

        // restricting global matrix into local submatrixes 
        bool restrict(const GlobalMatrix & M, LocalMatrix & M_k) const override
        { 
            this->_decomposition->restrict(M, M_k);
            return true;
        }



    private:    
        std::shared_ptr<dolfin::Function> _u;                               // solution
        std::shared_ptr<const dolfin::Form> _E;                             // energy 
        std::shared_ptr<const dolfin::Form> _g;                             // gradient
        std::shared_ptr<const dolfin::Form> _H;                             // hessian
        std::vector<std::shared_ptr<const dolfin::DirichletBC>> _bcs;       // BC 

        mutable PetscInt local_dim; 
        std::shared_ptr<FenicsDecomposition> _decomposition;
        SizeType mpi_size = mpi_world_size();
        SizeType mpi_rank = mpi_world_rank();
    };

}




#endif  //FENICS_UTOPIA_GLOBAL_LOCAL_HPP
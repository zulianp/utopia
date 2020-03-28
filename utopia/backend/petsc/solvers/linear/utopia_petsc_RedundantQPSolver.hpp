#ifndef UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP
#define UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP


#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_MPRGP.hpp"

namespace utopia {


    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class RedundantQPSolver {
    public:
        RedundantQPSolver() {
            static_assert(Backend < HOMEMADE, "RedundantLinearSolver not implemented for this backend");
        }

    };


    template<typename Matrix, typename Vector>
    class RedundantQPSolver<Matrix, Vector, PETSC> final : public OperatorBasedQPSolver<Matrix, Vector>
    {
        typedef UTOPIA_SCALAR(Vector)       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;
        using QPSolver         = utopia::OperatorBasedQPSolver<Matrix, Vector>;

        public:
            RedundantQPSolver(const std::shared_ptr<QPSolver> &qp_solver): 
            n_solves_(2), init_(false), qp_solver_(qp_solver)
            {

            }

            ~RedundantQPSolver()
            {
                VecScatterDestroy(&scatterout); 
                VecScatterDestroy(&scatterin); 
                PetscSubcommDestroy(&psubcomm);
            }


            RedundantQPSolver * clone() const override
            {
                return new RedundantQPSolver(*this);
            }


            void number_of_parallel_solves(const SizeType & number)
            {
                n_solves_ = number; 
            }

            void update(const Operator<Vector> &A) override
            {
                std::cout<<"------ not implemented \n"; 
            }            


            bool solve(const Operator<Vector> &A, const Vector &rhs, Vector &sol) override
            {
                std::cout<<"------- not implemented ---- \n"; 
                return false; 
            }


            bool solve(const Matrix &A, const Vector &rhs, Vector &sol) override
            {
                if(init_==false){
                    init_redundant(A, rhs); 
                }
                
                global_to_sub(sol, sol_sub, rhs, rhs_sub);

            
                // // // // // // // // // // // // //  lets assume, there will be solve here... // // // // // // // // 
                auto QP_solver = std::make_shared<utopia::MPGRP<Matrix, Vector> >();
                
                Vector lb = -9e9*sol_sub; 
                Vector ub = 9e9*sol_sub; 
                auto box = make_box_constaints(make_ref(lb), make_ref(ub)); 
                QP_solver->max_it(10);
                QP_solver->verbose(true); 

                QP_solver->init_memory(sol_sub.comm(), size(sol_sub).get(0), local_size(sol_sub).get(0)); 
                QP_solver->aux_solve(pmats, rhs_sub, sol_sub, box);

                

                
                // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
                // std::cout<<"---------------- \n"; 
                sub_to_global(sol_sub, sol); 

                


                return false; 
            }



        void init_redundant(const Matrix A, const Vector & x)
        {
            auto comm = A.comm().get(); 

            PetscSubcommCreate(comm, &psubcomm);
            PetscSubcommSetNumber(psubcomm, n_solves_);
            PetscSubcommSetType(psubcomm, PETSC_SUBCOMM_CONTIGUOUS);       

            MPI_Comm subcomm = PetscSubcommChild(psubcomm);

            pmats.destroy(); 
            MatCreateRedundantMatrix(raw_type(A), psubcomm->n, subcomm, MAT_INITIAL_MATRIX, &raw_type(pmats));

            // construct scatters 
            PetscInt       mstart, mend, mlocal, M, mloc_sub;

            /* get working vectors sol_sub and rhs_sub */
            MatCreateVecs(raw_type(pmats),  &raw_type(sol_sub), & raw_type(rhs_sub));


            MatGetLocalSize(raw_type(pmats), &mloc_sub, NULL);
            VecCreateMPI(PetscSubcommContiguousParent(psubcomm), mloc_sub, PETSC_DECIDE, & raw_type(sol_dup));
            VecCreateMPI(PetscSubcommContiguousParent(psubcomm), mloc_sub, PETSC_DECIDE, & raw_type(rhs_dup));


            // create scatters 
            IS       is1, is2;
            PetscInt *idx1,*idx2,i,j,k;

            VecGetSize(raw_type(x), &M);
            VecGetOwnershipRange(raw_type(x), &mstart, &mend);
            mlocal = mend - mstart;
            PetscMalloc2(psubcomm->n*mlocal, &idx1, psubcomm->n*mlocal, &idx2);
            j    = 0;

            for (k=0; k<psubcomm->n; k++) {
              for (i=mstart; i<mend; i++) {
                idx1[j]   = i;
                idx2[j++] = i + M*k;
              }
            }

            ISCreateGeneral(comm, psubcomm->n*mlocal, idx1, PETSC_COPY_VALUES, &is1);
            ISCreateGeneral(comm, psubcomm->n*mlocal, idx2, PETSC_COPY_VALUES, &is2);
            VecScatterCreate(raw_type(x), is1, raw_type(sol_dup), is2, &scatterin);
            ISDestroy(&is1);
            ISDestroy(&is2);


            ISCreateStride(comm, mlocal, mstart+ psubcomm->color*M, 1, &is1);
            ISCreateStride(comm, mlocal, mstart, 1, &is2);
            VecScatterCreate(raw_type(sol_dup), is1, raw_type(x), is2, &scatterout);
            ISDestroy(&is1);
            ISDestroy(&is2);
            PetscFree2(idx1, idx2);



            init_ = true; 
        }


        void global_to_sub(const Vector & sol, Vector & sol_sub, const Vector & rhs, Vector & rhs_sub)
        {
            // scatter solution, rhs 
            PetscScalar    *array_sol, *array_rhs;

            VecScatterBegin(scatterin, raw_type(sol), raw_type(sol_dup), INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatterin, raw_type(sol), raw_type(sol_dup), INSERT_VALUES, SCATTER_FORWARD);

            VecGetArray(raw_type(sol_dup), &array_sol);
            VecPlaceArray(raw_type(sol_sub),  (const PetscScalar*)array_sol);
            
            
            VecScatterBegin(scatterin, raw_type(rhs), raw_type(rhs_dup), INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatterin, raw_type(rhs), raw_type(rhs_dup), INSERT_VALUES, SCATTER_FORWARD);

            
            VecGetArray(raw_type(rhs_dup), &array_rhs);
            VecPlaceArray(raw_type(rhs_sub),  (const PetscScalar*)array_rhs);
            
        }


        void sub_to_global(const Vector & sol_sub, Vector & sol)
        {
            PetscScalar    *array_sol; 

            /* place ysub's local array into ydup */
            VecGetArray(raw_type(sol_sub), &array_sol);
            VecPlaceArray(raw_type(sol_dup), (const PetscScalar*)array_sol);


            /* scatter ydup to y */
            VecScatterBegin(scatterout, raw_type(sol_dup), raw_type(sol), INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatterout, raw_type(sol_dup), raw_type(sol), INSERT_VALUES, SCATTER_FORWARD);


            VecResetArray(raw_type(sol_dup));
            VecRestoreArray(raw_type(sol_sub), &array_sol);
        }



        private:
            SizeType n_solves_;
            bool init_; 
            std::shared_ptr<QPSolver> qp_solver_;     /*!< Linear solver parameters. */        


            PetscSubcomm        psubcomm;    
            Matrix              pmats;  // check for update 

            Vector         sol_sub, rhs_sub, sol_dup, rhs_dup;
            VecScatter     scatterin, scatterout; 

    };
}

#endif //UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP
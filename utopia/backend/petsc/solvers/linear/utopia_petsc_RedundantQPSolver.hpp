#ifndef UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP
#define UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP


#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"

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
            RedundantQPSolver(const std::shared_ptr<QPSolver> &qp_solver)
            : n_solves_(2), qp_solver_(qp_solver)
            {

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
                // std::cout<<"Size: "<< size(A) << "  \n"; 
                // std::cout<<"Size: "<< size(rhs) << "  \n"; 
                // std::cout<<"Size: "<< size(sol) << "  \n"; 


                // TODO:: check if we are in parallel, otherwise just run QP solver 

                // creating sub-comunicators 
                MPI_Comm            comm,subcomm;
                PetscSubcomm        psubcomm;


                comm = A.comm().get(); 
                
                PetscSubcommCreate(A.comm().get(), &psubcomm);
                PetscSubcommSetNumber(psubcomm, n_solves_);
                PetscSubcommSetType(psubcomm, PETSC_SUBCOMM_CONTIGUOUS);       

                subcomm = PetscSubcommChild(psubcomm);

                // PetscInt size; 
                // MPI_Comm_size(subcomm, &size);
                // PetscPrintf(subcomm, "size %d\n", (int)size);

                // destroy P mat 
                // construct redundant matrix 
                // reuse them for RMTR coarse grid 
                Matrix pmats; 
                MatCreateRedundantMatrix(raw_type(A), psubcomm->n, subcomm, MAT_INITIAL_MATRIX, &raw_type(pmats));
                // disp(pmats); 

            
                // construct scatters 
                PetscInt       mstart, mend, mlocal, M, mloc_sub;
                Vector         sol_sub, rhs_sub;            /* vectors of a subcommunicator to hold parallel vectors of PetscObjectComm((PetscObject)pc) */
                Vec            xdup,ydup;            /* parallel vector that congregates sol_sub or rhs_sub facilitating vector scattering */
                VecScatter     scatterin,scatterout; /* scatter used to move all values to each processor group (subcommunicator) */                


                /* get working vectors sol_sub and rhs_sub */
                MatCreateVecs(raw_type(pmats),  &raw_type(sol_sub), & raw_type(rhs_sub));


                /* create working vectors xdup and ydup.
                 xdup concatenates all xsub's contigously to form a mpi vector over dupcomm  (see PetscSubcommCreate_interlaced())
                 ydup concatenates all rhs_sub and has empty local arrays because rhs_sub's arrays will be place into it.
                 Note: we use communicator dupcomm, not PetscObjectComm((PetscObject)pc)! */
                MatGetLocalSize(raw_type(pmats), &mloc_sub, NULL);
                VecCreateMPI(PetscSubcommContiguousParent(psubcomm), mloc_sub, PETSC_DECIDE, &xdup);
                VecCreateMPI(PetscSubcommContiguousParent(psubcomm), mloc_sub, PETSC_DECIDE, &ydup);
                // VecCreateMPIWithArray(PetscSubcommContiguousParent(psubcomm), 1, mloc_sub, PETSC_DECIDE, NULL, &ydup);                





            /* create vecscatters */
            // if (!scatterin) 
            { 
                IS       is1,is2;
                PetscInt *idx1,*idx2,i,j,k;
                Vec x; 

                MatCreateVecs(raw_type(A), &x, 0);
                VecGetSize(x, &M);
                VecGetOwnershipRange(x, &mstart, &mend);
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
                VecScatterCreate(x, is1, xdup, is2, &scatterin);
                ISDestroy(&is1);
                ISDestroy(&is2);

                /* Impl below is good for PETSC_SUBCOMM_INTERLACED (no inter-process communication) and PETSC_SUBCOMM_CONTIGUOUS (communication within subcomm) */
                ISCreateStride(comm, mlocal, mstart+ psubcomm->color*M, 1, &is1);
                ISCreateStride(comm, mlocal, mstart, 1, &is2);
                VecScatterCreate(xdup, is1, x, is2, &scatterout);
                ISDestroy(&is1);
                ISDestroy(&is2);
                PetscFree2(idx1, idx2);
                VecDestroy(&x);
            }

            // scatter solution, rhs 
            PetscScalar    *array_sol, *array_rhs;

            /* scatter x to xdup */
            VecScatterBegin(scatterin, raw_type(sol), xdup, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatterin, raw_type(sol), xdup, INSERT_VALUES, SCATTER_FORWARD);

            /* place xdup's local array into xsub */
            VecGetArray(xdup, &array_sol);
            VecPlaceArray(raw_type(sol_sub),  (const PetscScalar*)array_sol);
            // disp(sol_sub); 


            /* scatter x to xdup */
            VecScatterBegin(scatterin, raw_type(rhs), ydup, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatterin, raw_type(rhs), ydup, INSERT_VALUES, SCATTER_FORWARD);


            /* place xdup's local array into xsub */
            VecGetArray(ydup, &array_rhs);
            VecPlaceArray(raw_type(rhs_sub),  (const PetscScalar*)array_rhs);
            // disp(rhs_sub);             


            // qp_solver_->solve(pmats, rhs_sub, sol_sub); 
            // std::cout<<"1Size: "<< size(pmats) << "  \n"; 
            // std::cout<<"1Size: "<< size(rhs_sub) << "  \n"; 
            // std::cout<<"1Size: "<< size(sol_sub) << "  \n"; 

            // // // // // // // // // // //  lets assume, there will be solve here... // // // // // // // // 
            Vector bla = 0.0* sol_sub; 
            bla = pmats * rhs_sub;
            disp(bla); 
            // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 

            
            /* place ysub's local array into ydup */
            VecGetArray(raw_type(sol_sub), &array_sol);
            VecPlaceArray(xdup, (const PetscScalar*)array_sol);


            /* scatter ydup to y */
            VecScatterBegin(scatterout, xdup, raw_type(sol), INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(scatterout, xdup, raw_type(sol), INSERT_VALUES, SCATTER_FORWARD);


            VecResetArray(xdup);
            VecRestoreArray(raw_type(sol_sub), &array_sol);



            disp(sol, "sol"); 
            // TODO make sure that everything is deleted 



            std::cout<<"----- here ---- \n";
            exit(0);

            return false; 
        }






        private:
            SizeType n_solves_;
            std::shared_ptr<QPSolver> qp_solver_;     /*!< Linear solver parameters. */            
    };
}

#endif //UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP
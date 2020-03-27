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

                // creating sub-comunicators 
                MPI_Comm            comm,subcomm;
                PetscSubcomm        psubcomm;
                
                PetscSubcommCreate(A.comm().get(), &psubcomm);
                PetscSubcommSetNumber(psubcomm, n_solves_);
                PetscSubcommSetType(psubcomm, PETSC_SUBCOMM_CONTIGUOUS);       

                subcomm = PetscSubcommChild(psubcomm);

                // PetscInt size; 
                // MPI_Comm_size(subcomm, &size);
                // PetscPrintf(subcomm, "size %d\n", (int)size);

                // destroy P mat 
                // construct redundant matrix 
                Matrix pmats; 
                MatCreateRedundantMatrix(raw_type(A), psubcomm->n, subcomm, MAT_INITIAL_MATRIX, &raw_type(pmats));
                disp(pmats); 

                // MatCreateVecs(red->pmats,&red->xsub,&red->ysub);





                // std::cout<<"----- here ---- \n";
                exit(0);

                return false; 
            }






        private:
            SizeType n_solves_;
            std::shared_ptr<QPSolver> qp_solver_;     /*!< Linear solver parameters. */            
    };
}

#endif //UTOPIA_PETSC_REDUNDANT_QP_SOLVER_HPP
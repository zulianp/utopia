#ifndef UTOPIA_PETSC_REDUNDANT_LINEAR_SOLVER_HPP
#define UTOPIA_PETSC_REDUNDANT_LINEAR_SOLVER_HPP


#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"

namespace utopia {
    template<typename Matrix, typename Vector>
    class RedundantLinearSolver<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> 
    {
        typedef UTOPIA_SCALAR(Vector)       Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector)    SizeType;

        public:
            RedundantLinearSolver(const std::string &sub_preconditioner = "lu")
            : KSPSolver<Matrix, Vector, PETSC>()
            {
                KSPSolver<Matrix, Vector, PETSC>::pc_type("redundant");
                KSPSolver<Matrix, Vector, PETSC>::ksp_type("preonly");

                this->pc_type(sub_preconditioner); 
            }

            RedundantLinearSolver * clone() const override
            {
                return new RedundantLinearSolver(*this);
            }


            void number_of_parallel_solves(const SizeType & number)
            {
                PC pc_redundant; 
                KSPGetPC(this->implementation(), &pc_redundant);                    
                PCRedundantSetNumber(pc_redundant, number); 
            }


            void pc_type(const std::string &pc_type) override
            {
                PC pc_redundant; 
                KSPGetPC(this->implementation(), &pc_redundant);      

                KSP innerksp; PC inner_pc; 
                PCRedundantGetKSP(pc_redundant, &innerksp); 
                KSPGetPC(innerksp, &inner_pc);                
                PCSetType(inner_pc, pc_type.c_str()); 
            }

            void ksp_type(const std::string &ksp_type) override
            {
                PC pc_redundant; 
                KSPGetPC(this->implementation(), &pc_redundant);      

                KSP innerksp;
                PCRedundantGetKSP(pc_redundant, &innerksp); 
                KSPSetType(innerksp, ksp_type.c_str()); 
            }
    };
}

#endif //UTOPIA_PETSC_REDUNDANT_LINEAR_SOLVER_HPP
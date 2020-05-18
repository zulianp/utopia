#ifndef UTOPIA_PETSC_GMRES_HPP
#define UTOPIA_PETSC_GMRES_HPP

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_petsc_KSPSolver.hpp"

namespace utopia {
    template <typename Matrix, typename Vector>
    class GMRES<Matrix, Vector, PETSC> final : public KSPSolver<Matrix, Vector, PETSC> {
    public:
        GMRES(const std::string &preconditioner = "jacobi") : KSPSolver<Matrix, Vector, PETSC>() {
            this->pc_type(preconditioner);
            this->ksp_type("gmres");
        }

        GMRES *clone() const override { return new GMRES(*this); }
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_GMRES_HPP

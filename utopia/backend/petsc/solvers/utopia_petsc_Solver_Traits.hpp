#ifndef UTOPIA_PETSC_SOLVER_TRAITS_HPP
#define UTOPIA_PETSC_SOLVER_TRAITS_HPP

#include "utopia_Traits.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"

namespace utopia {

    template <typename Matrix, typename Vector, int Backend>
    class Traits<KSPSolver<Matrix, Vector, Backend> > : public Traits<Matrix> {};
}  // namespace utopia

#endif  // UTOPIA_PETSC_SOLVER_TRAITS_HPP

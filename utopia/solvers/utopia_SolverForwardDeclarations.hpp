#ifndef UTOPIA_SOLVERFORWARDDECLARATIONS_HPP
#define UTOPIA_SOLVERFORWARDDECLARATIONS_HPP

namespace utopia {

    template<class Matrix, class Vector>
    class IterativeSolver;

    template<class Matrix, class Vector>
    class LinearSolver;

    template<typename Matrix, typename Vector>
    class DirectSolver;

    template<class Vector>
    class MatrixFreeLinearSolver;

    template<class Matrix, class Vector>
    class OperatorBasedLinearSolver;

    template<class Vector>
    class Preconditioner;

    template<class Matrix, class Vector>
    class PreconditionedSolver;

}

#endif //UTOPIA_SOLVERFORWARDDECLARATIONS_HPP

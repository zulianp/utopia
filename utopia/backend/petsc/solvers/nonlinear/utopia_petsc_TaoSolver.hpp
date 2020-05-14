#ifndef UTOPIA_PETSC_TAO_SOLVER_HPP
#define UTOPIA_PETSC_TAO_SOLVER_HPP

#include "utopia_BoxConstraints.hpp"
#include "utopia_NonLinearSolver.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_KSPSolver.hpp"
#include "utopia_petsc_Types.hpp"
#include "utopia_Function.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"
#include "utopia_petsc_build_ksp.hpp"
#include "utopia_NonLinearSmoother.hpp"

#include <string>

namespace utopia {

    template<class Matrix, class Vector>
    class TaoSolver final : public NewtonBase<Matrix, Vector>,
                            // public NonLinearSmoother<Matrix, Vector>, //Maybe removing rhs from interface??
                            public VariableBoundSolverInterface<Vector>,
                            public virtual Clonable {
    public:
        TaoSolver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &linear_solver);
        TaoSolver();
        ~TaoSolver() override;

        void set_type(const std::string &type);
        void read(Input &in) override;
        void print_usage(std::ostream &os) const override;

        bool solve(Function<Matrix, Vector> &fun, Vector &x) override;
        bool smooth(Function<Matrix, Vector> &fun, Vector &x);// override;
        TaoSolver * clone() const override;


    private:
        class Impl;
        std::unique_ptr<Impl> impl_;

        void init(Function<Matrix, Vector> &fun, Vector &x);
        void set_function(Function<Matrix, Vector> &fun);
        bool get_ksp(KSP *ksp);
    };
}

#endif //UTOPIA_PETSC_TAO_SOLVER_HPP

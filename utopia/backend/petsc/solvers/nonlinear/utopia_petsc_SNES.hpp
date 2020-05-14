#ifndef UTOPIA_SNES_SOLVER_HPP
#define UTOPIA_SNES_SOLVER_HPP

#include "utopia_Core.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_Function.hpp"
#include "utopia_NonLinearSolver.hpp"
// #include "utopia_petsc.hpp"
#include "utopia_petsc_build_ksp.hpp"
#include "utopia_VariableBoundSolverInterface.hpp"

#include <algorithm>
#include <petscpc.h>
#include <petscksp.h>
#include <petscsys.h>
#include <petscsnes.h>

namespace utopia {
    // this should be potentially nonlinear solver ....
    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class SNESSolver {};

    template<typename Matrix, typename Vector>
    class SNESSolver<Matrix, Vector, PETSC> : public NewtonBase<Matrix, Vector>,
                                              public NonLinearSmoother<Matrix, Vector>,
                                              public virtual VariableBoundSolverInterface<Vector>,
                                              public virtual Clonable //FIXME all non-linear solvers should be clonable
{
    public:
        using Scalar   = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout   = typename Traits<Vector>::Layout;

        typedef typename NewtonBase<Matrix, Vector>::Solver  LinearSolver;
        typedef utopia::NonLinearSmoother<Matrix, Vector>    Smoother;
        typedef utopia::NewtonBase<Matrix, Vector>           NonLinearSolver;
        typedef utopia::Function<Matrix, Vector>             Function;

        SNESSolver(const std::shared_ptr<LinearSolver> &linear_solver = std::shared_ptr<LinearSolver>(),
                   const std::vector<std::string> snes_types    = {"newtonls", "newtontr", "nrichardson", "ksponly", "vinewtonrsls", "vinewtonssls", "ngmres", "qn", "shell", "ngs", "ncg", "fas", "ms", "anderson"});

        SNESSolver *clone() const override;

        ~SNESSolver() override;

        void read(Input &in) override;

        void print_usage(std::ostream &os) const override;

        virtual void set_snes_type(const std::string & type);

        bool solve(Function &fun, Vector &x) override;

        bool smooth(Function &fun, Vector &x, const Vector &rhs) override;

        void set_line_search_type(const SNESLineSearchType ls_type);

    protected:

        std::string SNES_type_;                                  /*!< Choice of snes types. */
        const std::vector<std::string> SNES_types;              /*!< Valid options for SNES solver types. */
        SNESLineSearchType line_search_type_;

        virtual void set_snes_options(SNES & snes,
                                      const Scalar & atol,
                                      const Scalar & rtol,
                                      const Scalar & stol,
                                      const SizeType & max_it);


        virtual void set_ksp(SNES & snes);

    private:

        void set_variable_bounds(SNES &snes, const Layout &layout);

        // this function is expensive - wrt to convert
        // however, in utopia is not ready for pure wrap implementation
        void setup_assembly_routines(SNES & snes, Function & fun, const Vector &x);
    };
}

#endif // UTOPIA_SNES_SOLVER_HPP

#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Base.hpp"
#include "utopia_Testing.hpp"

#include "utopia_InputParameters.hpp"

#include "utopia_Bratu1D.hpp"
#include "utopia_SubCommUnitTest.hpp"

#include "utopia_PseudoTimeStepper.hpp"
#include "utopia_TwoFieldAlternateMinimization.hpp"
#include "utopia_TwoFieldSPIN.hpp"

#include <memory>

#include "utopia_Poisson1D.hpp"

#ifdef UTOPIA_WITH_TRILINOS
#include "utopia_Tpetra_Operator.hpp"
#include "utopia_trilinos_GMRES.hpp"
#endif

#ifdef UTOPIA_WITH_PETSC
#include "utopia_petsc_GMRES.hpp"
#endif

using namespace utopia;

template <class Matrix, class Vector = typename Traits<Matrix>::Vector>
class TimeDependentPoisson1D : public TimeDependentFunction<Matrix, Vector> {
public:
    using Traits = utopia::Traits<Vector>;
    using Scalar = typename Traits::Scalar;
    using SizeType = typename Traits::SizeType;

    TimeDependentPoisson1D(const int n) : fun_(std::make_shared<Poisson1D<Matrix>>(n)) {}

    bool hessian(const Vector &x, Matrix &H) const override { return fun_->hessian(x, H); }
    bool value(const Vector &x, Scalar &value) const override { return fun_->value(x, value); }
    bool gradient(const Vector &x, Vector &g) const override { return fun_->gradient(x, g); }
    void create_vector(Vector &x) const override { fun_->create_vector(x); }
    bool update(const Scalar t) override { return true; }

    std::shared_ptr<Function<Matrix, Vector>> fun_;
};

template <class Matrix, class Vector>
class FieldSplitTestTest : public SubCommUnitTest<Vector> {
public:
    using Traits = utopia::Traits<Vector>;
    using Scalar = typename Traits::Scalar;
    using SizeType = typename Traits::SizeType;
    using Comm = typename Traits::Communicator;

    void test_aternate_minimization() {
        int n = 10;
        bool verbose = false;

        // Split-fields
        auto f1 = std::make_shared<Poisson1D<Matrix>>(n);
        auto f2 = std::make_shared<Poisson1D<Matrix>>(n);

        // Global optimization problem (coupled fields)
        auto c12 = std::make_shared<Poisson1D<Matrix>>(n);

        auto nls1 = std::make_shared<Newton<Matrix>>();
        auto nls2 = std::make_shared<Newton<Matrix>>();

        nls1->verbose(verbose);
        nls2->verbose(verbose);

        TwoFieldAlternateMinimization<Matrix> tfa(nls1, nls2);

        tfa.set_field_functions(f1, f2);

        auto f1_to_c12 = [](const Vector &in, Vector &out) { out = in; };
        auto f2_to_c12 = f1_to_c12;
        auto c12_to_f1 = f1_to_c12;
        auto c12_to_f2 = f1_to_c12;

        tfa.set_transfers(f1_to_c12, f2_to_c12, c12_to_f1, c12_to_f2);

        tfa.verbose(verbose);

        Vector x = c12->initial_guess();
        tfa.solve(*c12, x);
    }

    auto cg() -> std::shared_ptr<ConjugateGradient<Matrix, Vector, HOMEMADE>> {
        auto ret = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
        ret->apply_gradient_descent_step(true);
        return ret;
    }

    auto mf_solver() -> std::shared_ptr<MatrixFreeLinearSolver<Vector>> {
        if constexpr (utopia::Traits<Matrix>::Backend == TRILINOS) {
            auto ret = std::make_shared<GMRES<Matrix, Vector>>();
            ret->verbose(true);
            return ret;
        } else if constexpr (utopia::Traits<Matrix>::Backend == PETSC) {
            auto ret = std::make_shared<utopia::BiCGStab<Matrix, Vector, HOMEMADE>>();
            return ret;
        }
    }

    void test_SPIN() {
        int n = 10;
        bool verbose = true;

        // Split-fields
        auto f1 = std::make_shared<Poisson1D<Matrix>>(n);
        auto f2 = std::make_shared<Poisson1D<Matrix>>(n);

        // Global optimization problem (coupled fields)
        auto c12 = std::make_shared<Poisson1D<Matrix>>(n);

        auto nls1 = std::make_shared<Newton<Matrix>>(cg());
        auto nls2 = std::make_shared<Newton<Matrix>>(cg());

        // Global linear solver
        auto lsc12 = mf_solver();
        // auto lsc12 = cg();

        nls1->verbose(verbose);
        nls2->verbose(verbose);

        TwoFieldSPIN<Matrix> tfa(lsc12, nls1, nls2);
        tfa.set_field_functions(f1, f2);
        // tfa.additive_precond(false);

        auto f1_to_c12 = [](const Vector &in, Vector &out) { out = in; };
        auto f2_to_c12 = f1_to_c12;
        auto c12_to_f1 = f1_to_c12;
        auto c12_to_f2 = f1_to_c12;

        tfa.set_transfers(f1_to_c12, f2_to_c12, c12_to_f1, c12_to_f2);

        tfa.verbose(verbose);

        // if (verbose) {
        //     tfa.verbosity_level(VERBOSITY_LEVEL_DEBUG);
        // }

        Vector x = c12->initial_guess();
        tfa.solve(*c12, x);
    }

    void test_incremental_loading() {
        bool verbose = false;

        auto nlsolver = std::make_shared<Newton<Matrix>>();
        nlsolver->verbose(verbose);

        PseudoTimeStepper<Matrix, Vector> time_stepper(nlsolver);
        auto poisson = std::make_shared<TimeDependentPoisson1D<Matrix>>(10);

        Vector x;
        poisson->create_vector(x);
        time_stepper.set_end_time(10);
        time_stepper.solve(*poisson, x);
        disp(x);
    }

    void run() {
        UTOPIA_RUN_TEST(test_aternate_minimization);
        UTOPIA_RUN_TEST(test_SPIN);
        UTOPIA_RUN_TEST(test_incremental_loading);
    }
};

void field_split() {
#ifdef UTOPIA_WITH_PETSC
    FieldSplitTestTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC

#ifdef UTOPIA_WITH_TRILINOS
    FieldSplitTestTest<TpetraMatrix, TpetraVector>().run();
#endif
}

UTOPIA_REGISTER_TEST_FUNCTION(field_split);

#include "utopia.hpp"
#include "utopia_Assert.hpp"
#include "utopia_Testing.hpp"

#include "utopia_InputParameters.hpp"

#include "utopia_Bratu1D.hpp"
#include "utopia_SubCommUnitTest.hpp"

#include "utopia_TwoFieldAlternateMinimization.hpp"
#include "utopia_TwoFieldSPIN.hpp"

#include <memory>

#include "utopia_Poisson1D.hpp"

using namespace utopia;

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

    void test_SPIN() {
        int n = 10;
        bool verbose = true;

        // Split-fields
        auto f1 = std::make_shared<Poisson1D<Matrix>>(n);
        auto f2 = std::make_shared<Poisson1D<Matrix>>(n);

        // Global optimization problem (coupled fields)
        auto c12 = std::make_shared<Poisson1D<Matrix>>(n);

        auto nls1 = std::make_shared<Newton<Matrix>>();
        auto nls2 = std::make_shared<Newton<Matrix>>();

        // Global linear solver
        auto lsc12 = std::make_shared<ConjugateGradient<Matrix, Vector, HOMEMADE>>();
        lsc12->apply_gradient_descent_step(true);

        nls1->verbose(verbose);
        nls2->verbose(verbose);

        TwoFieldSPIN<Matrix> tfa(lsc12, nls1, nls2);
        tfa.set_field_functions(f1, f2);
        tfa.additive_precond(false);

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

    void run() {
        UTOPIA_RUN_TEST(test_aternate_minimization);
        UTOPIA_RUN_TEST(test_SPIN);
    }
};

void field_split() {
#ifdef UTOPIA_WITH_PETSC
    FieldSplitTestTest<PetscMatrix, PetscVector>().run();
#endif  // UTOPIA_WITH_PETSC
}

UTOPIA_REGISTER_TEST_FUNCTION(field_split);

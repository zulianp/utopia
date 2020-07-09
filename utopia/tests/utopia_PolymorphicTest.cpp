#include <cassert>
#include <string>
#include "utopia.hpp"
#include "utopia_Testing.hpp"

#include "utopia_AbstractLinearSolver.hpp"
#include "utopia_AbstractVector.hpp"
#include "utopia_ObjectFactory.hpp"
#include "utopia_make_unique.hpp"

#ifdef WITH_PETSC
#include "utopia_petsc.hpp"
#include "utopia_petsc_impl.hpp"
#endif  // WITH_PETSC

namespace utopia {

// Register types
#ifdef WITH_PETSC

    // FIXME make it so that also includes utopia front-end solvers
    using PetscLinearSolver = utopia::KSPSolver<PetscMatrix, PetscVector>;

    UTOPIA_FACTORY_REGISTER_VECTOR(PetscVector);
    UTOPIA_FACTORY_REGISTER_MATRIX(PetscMatrix);
    UTOPIA_FACTORY_REGISTER_LINEAR_SOLVER(PetscLinearSolver);

#ifdef WITH_TRILINOS
    UTOPIA_FACTORY_REGISTER_VECTOR(TpetraVector);
#endif

    class PolymorphicTest final {
    public:
        using Scalar = typename Traits<PetscVector>::Scalar;
        using SizeType = typename Traits<PetscVector>::SizeType;

        using DefaultFactory = utopia::AlgebraFactory<Scalar, SizeType>;

        // base classes
        using AbstractVector = utopia::AbstractVector<Scalar, SizeType>;
        using AbstractMatrix = utopia::AbstractMatrix<Scalar, SizeType>;
        using AbstractLinearSolver = utopia::AbstractLinearSolver<Scalar, SizeType>;

        void convenience_wrapper() {
            const SizeType n = 10;

#ifdef WITH_TRILINOS
#ifdef UTOPIA_TPETRA_SIZE_TYPE
            // if types are the same TirlinosFactory == DefaultFactory
            InputParameters params;
            params.set("default-backend", "trilinos");
            DefaultFactory::init(params);

#endif  // UTOPIA_TPETRA_SIZE_TYPE
#endif  // WITH_TRILINOS

            auto x = DefaultFactory::new_vector();
            x->values(serial_layout(n), 2.0);

            auto m = unique_to_shared(DefaultFactory::new_matrix());

            if (m) {
                m->identity(serial_layout(n, n), 2.0);

                auto y = DefaultFactory::new_vector();
                m->apply(*x, *y);

                Scalar y_n = y->norm2();
                utopia_test_assert(approxeq(std::sqrt(n * 16.0), y_n, 1e-8));

                auto s = DefaultFactory::new_linear_solver();

                // InputParameters sol_params;
                // sol_params.set("verbose", true);
                // s->read(sol_params);

                // set initial guess to 0
                x->set(0.0);

                // or call s->solve(*m, *x, *y); if applied only once for m
                s->update(m);
                s->apply(*y, *x);

                auto r = DefaultFactory::new_vector();

                m->apply(*x, *r);
                r->axpy(-1.0, *y);

                Scalar norm_r = r->norm2();

                disp(norm_r);
                utopia_test_assert(norm_r < 1e-6);
            }
        }

        void run() { UTOPIA_RUN_TEST(convenience_wrapper); }
    };

    static void polymorphic() { PolymorphicTest().run(); }

    UTOPIA_REGISTER_TEST_FUNCTION(polymorphic);

#endif  // WITH_PETSC
}  // namespace utopia

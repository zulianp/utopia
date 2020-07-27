#include "utopia.hpp"

#include <algorithm>
#include <cmath>

using namespace utopia;

template <class Vector>
static void run(const int n) {
    // using Scalar = typename Traits<Vector>::Scalar;

    utopia::out() << "Each Utopia" << std::endl;
    Vector v;

    v.values(serial_layout(n), 1.);
    utopia::out() << backend_info(v).get_name() << std::endl;

    Chrono c;
    c.start();

    {
        auto v_view = local_view_device(v);
        parallel_for(
            local_range_device(v), UTOPIA_LAMBDA(const SizeType &i) { v_view.set(i, 10.); });
    }

    c.stop();
    utopia::out() << c << std::endl;
}

static void run_raw(const int n) {
    utopia::out() << "Raw" << std::endl;
    std::vector<double> v(n, 1.);

    Chrono c;
    c.start();
    for (auto &vi : v) {
        vi = 10.;
    }
    c.stop();
    utopia::out() << c << std::endl;
}

static void run_array(const int n) {
    utopia::out() << "Array" << std::endl;
    std::vector<double> v(n, 1.);

    auto a = &v[0];

    Chrono c;
    c.start();

    for (int i = 0; i < n; ++i) {
        a[i] = 10.;
    }

    c.stop();
    utopia::out() << c << std::endl;
}

static void run_stl_each(const int n) {
    utopia::out() << "ForEach" << std::endl;
    std::vector<double> v(n, 1.);

    Chrono c;
    c.start();

    std::for_each(std::begin(v), std::end(v), [](double &value) { value = 10.; });

    c.stop();
    utopia::out() << c << std::endl;
}

static void run_for(const int n) {
    utopia::out() << "For" << std::endl;
    std::vector<double> v(n, 1.);

    Chrono c;
    c.start();

    For<64>::apply(0, n, [&v](const std::size_t i) { v[i] = 10.; });

    c.stop();
    utopia::out() << c << std::endl;
}

template <class Vector>
static void run_access(const int n) {
    utopia::out() << "Utopia Access" << std::endl;
    Vector v;
    v.values(serial_layout(n), 1.);

    utopia::out() << backend_info(v).get_name() << std::endl;

    Chrono c;
    c.start();

    {
        Write<Vector> w_(v);
        auto r = range(v);
        For<>::apply(r.begin(), r.end(), [&v](const SizeType i) { v.set(i, 10.); });
    }

    c.stop();
    utopia::out() << c << std::endl;
}

#ifdef WITH_BLAS
static void run_access_blas(const int n) {
    utopia::out() << "Utopia Access Blas" << std::endl;
    BlasVectord v;
    v.values(serial_layout(n), 1.0);
    utopia::out() << backend_info(v).get_name() << std::endl;

    Chrono c;
    c.start();

    {
        Write<BlasVectord> w_(v);
        auto r = range(v);
        For<>::apply(r.begin(), r.end(), [&v](const SizeType i) { v.set(i, 10); });
    }

    c.stop();
    utopia::out() << c << std::endl;
}
#endif  // WITH_BLAS

#ifdef WITH_PETSC
static void run_access_petsc(const int n) {
    utopia::out() << "Utopia Access Petsc" << std::endl;
    PetscVector v;
    v.values(serial_layout(n), 1.0);
    utopia::out() << backend_info(v).get_name() << std::endl;

    Chrono c;
    c.start();

    {
        PetscScalar val = 10.;
        Write<PetscVector> w_(v);
        auto r = range(v);
        For<>::apply(
            r.begin(), r.end(), [&v, val](const PetscInt i) { VecSetValues(raw_type(v), 1, &i, &val, INSERT_VALUES); });
    }

    c.stop();
    utopia::out() << c << std::endl;
}
#endif  // WITH_PETSC

static void run_all(const int n) {
    run_raw(n);
    run_array(n);
    run_stl_each(n);
    run_for(n);
    // if it has compiled with blas or petsc WITH_BLAS or WITH_PETSC macros are available (if you want to make it
    // compile no matter the utopia installation)
#ifdef WITH_PETSC
    // run with petsc types
    run<PetscVector>(n);
    run_access<PetscVector>(n);
    run_access_petsc(n);
#endif  // WITH_PETSC

#ifdef WITH_BLAS
    // run with blas types
    run<BlasVectord>(n);
    run_access<BlasVectord>(n);
    run_access_blas(n);
#endif  // WITH_BLAS

#ifdef WITH_TRILINOS
    // run with trilinos types
    run<TpetraVectord>(n);
    run_access<TpetraVectord>(n);
#endif  // WITH_TRILINOS
}

// Run it with `./examples/example_loops -n 100000
int main(int argc, char **argv) {
    using namespace utopia;

    Utopia::Init(argc, argv);

    InputParameters params;
    params.init(argc, argv);
    long n = 1000;
    params.get("n", n);
    run_all(n);

    return Utopia::Finalize();
}

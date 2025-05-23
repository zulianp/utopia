#include <iostream>
#include "utopia.hpp"

// Utopia may use different backend libraries
// to define a vector, a template class allows
// to write the code once for all vector types.
template <class Vector>
void sum_of_two_vectors() {
    using namespace utopia;

    // We define the type of the communicator. The communicator
    // allows different processes to communicate between them,
    // sharing local information for their computation.
    using Comm = typename Traits<Vector>::Communicator;

    // SyzeType is a type for indicies, in this case is used
    // for n_locals and n_global, which we will see in a second.
    using SizeType = typename Traits<Vector>::SizeType;

    // The element of a vector are of type Scalar.
    using Scalar = typename Traits<Vector>::Scalar;

    // We declare the vectors.
    Vector a, b;

    // Default initialisation of the communicator.
    Comm comm = Comm::get_default();

    // For each process we declare the number of local
    // entries that we want to use to compute our calculations.
    // In this case, our vector will be 10 x comm.size().
    SizeType n_local = 10;

    // We compute the global size of our vectors.
    SizeType n_global = n_local * comm.size();

    // A layout allows to describe the distribution of the data.
    // with n_local we define how many local entries we want
    // and with n_global how many global entries we have in total.
    auto l = layout(comm, n_local, n_global);

    // Initialise the vector `a` with uniform value `5`.
    a.values(l, 5);

    // Initialise the vector `b` with uniform value `2`.
    b.values(l, 2);

    // ----- Initialising vectors with different values ------

    // We have two way to initialise vectors, that is with local
    // or global indexing.
    // A local initialisation involves a fixed range for the entries. Each
    // entry may have a different processor.

    // You have a device (GPU, CPU) or a host (always CPU).
    // For using local indexing, we can do, for example:
    // auto a_view = local_view_device(a);.
    // Fo using global indexing, instead, with
    // auto a_view = view_device(a);.

    // You can run the programm to see the difference in terms
    // of output. The first output is with local indexing,
    // the second is with global indexing.

    // In both cases, we use the 'set(index, value)' function
    // to set a value at a certain index.
    // We can also use 'get(index)'' to get the value at certain
    // index.
    // the 'disp(vector/matrix)' function allows to visualise
    // the content of a vector or a matrix.

    {
        auto a_view = local_view_device(a);
        parallel_for(
            local_range_device(a), UTOPIA_LAMBDA(const SizeType& i) {
                const auto ai = a_view.get(i);
                a_view.set(i, i);
            });
    }

    // An expression object for sum of two vectors is created here, and ..
    auto c_expr = a + b;

    // evaluated here
    Vector c = c_expr;

    // Display the content.
    disp(c);

    {
        auto a_view = view_device(a);
        parallel_for(
            range_device(a), UTOPIA_LAMBDA(const SizeType& i) {
                const Scalar val = i;
                // a more complex example.
                // const Scalar val = (i == 0) ? 1e-14 : ((i < n / 2.0) ? -i : i);
                a_view.set(i, i);
            });
    }

    // We return the sum of the two vectors.
    auto d_expr = a + b;
    Vector d = d_expr;

    disp(d);

    // ---------------------------------------------------------------------
    // -------------- Wrong way to initialise a vector --------------------

    // This initialisation does not work since it is not performed
    // in parallel.
    // a_view.set(0, 1);
    // a_view.set(1, 10);
    // a_view.set(2, 33);

    // This does not work either since it is also not performed in
    // parallel.
    // auto j = 0;
    // for (auto i = 0; i < n_local; ++i) {
    //     a_view.set(i, j);
    //     j++;
    // }

    // ---------------------------------------------------------------------
}

// Use 'make -j4 complete' to compile.
// Then run with, for example, 4 processor with
// mpirun  -n 4 ./examples/your_example'
int main(int argc, char** argv) {
    using namespace utopia;

// Checking which backend is active.
#ifdef UTOPIA_ENABLE_PETSC
    using MatrixT = PetscMatrix;
    using VectorT = PetscVector;
#else
#ifdef UTOPIA_ENABLE_TRILINOS
    using MatrixT = TpetraMatrixd;
    using VectorT = TpetraVectord;
#else
    using MatrixT = BlasMatrixd;
    using VectorT = BlasVectord;
#endif
#endif

    Utopia::Init(argc, argv);

    sum_of_two_vectors<VectorT>();

    return Utopia::Finalize();
}

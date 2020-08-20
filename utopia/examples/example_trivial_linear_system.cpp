
#include <iostream>
#include <sstream>

#include "utopia.hpp"

// In this simple example we solve the following equation A x = b.
// A in this case is a simple Identity matrix = 1.
// So in the end we should have that x = b.

template <class Matrix, class Vector>
void example_trivial_linear_system(utopia::Input& in) {
    using namespace utopia;
    using Comm = typename Traits<Vector>::Communicator;
    using SizeType = typename Traits<Vector>::SizeType;

    // Instantiate the Matrix and vectors.
    Matrix A;
    Vector b, x;

    // Setup the communicator.
    Comm& comm = Comm::get_default();

    // Set default values.
    SizeType n_local = 10;
    auto value = 2;

    // Give different values to value and n_local from command line.
    in.get("n_local", n_local);
    in.get("value", value);

    // Setup the layout of the vectors and matrices.
    // Layouts allow us to distribute the data in according to the number of local entries n_local.
    // Note that the number of rows of the matrix and the size of the vectors are the same.
    SizeType n_global = n_local * comm.size();
    auto l = layout(comm, n_local, n_global);

    // Initialize the vector with the desired layout, and set all entries to zero.
    b.values(l, 0);

    {
        // Here we access the locally index view of b. Index goes from 0 to n_local-1.
        auto b_view = local_view_device(b);

        // For loop that iterates through the elements of the vector from index
        // 0 = r.begin() till index r.end() = n_local (excluded).
        // Change the entries of the vector at the local indexes 0 and 9.
        // This loop does not take advantage of parallelism and it is not hardware portable
        // (i.e., it will not work if the memory of b is stored on the GPU).
        // for (SizeType i = r.begin(); i != r.end(); ++i) {
        //     if (i > 0 && i < n_local - 1) {
        //         b.set(i, 0);
        //     }
        // }

        // Parallel for loop.
        parallel_for(
            local_range_device(b), UTOPIA_LAMBDA(const SizeType& i) {
                if (i == 0 || i == n_local - 1) {
                    b_view.set(i, value);
                }
            });
    }

    A.identity(square_matrix_layout(l), 1);

    // Example of assembling the identity matrix on the host.
    // {
    //     // Scoped lock for allowing to write data from the host.
    //     // When `w` is destroyed the memory is transferred to the device.
    //     // For parallel matrices the non-local data is communicated to the owner processes.
    //     Write<Matrix> w(A);
    //
    //     // r in a serial 3X3 matrix has range [0, 3), the row_range describes which local rows the matrix owns.
    //     auto r = row_range(A);
    //     for (SizeType i = r.begin(); i != r.end(); ++i) {
    //         A.set(i, i, 1);
    //     }
    // }

    // Solve the equation A x = b.
    solve(A, b, x);

    // Display the solution vector x.
    disp(x);
}

// Run it with `./examples/example_trivial_linear_system`
// or `./examples/example_trivial_linear_system -n_local 100 -value 23`
int main(const int argc, char* argv[]) {
    using namespace utopia;

    // Use certain Matrix or Vector depending
    // on what back-end we have installed on our system.
#ifdef UTOPIA_WITH_PETSC
    using MatrixT = PetscMatrix;
    using VectorT = PetscVector;
#else
#ifdef WITH_TRILINOS
    using MatrixT = TpetraMatrixd;
    using VectorT = TpetraVectord;
#else
    using MatrixT = BlasMatrixd;
    using VectorT = BlasVectord;
#endif
#endif

    // Initialize utopia...
    Utopia::Init(argc, argv);

    InputParameters params;
    params.init(argc, argv);

    // Run example
    example_trivial_linear_system<MatrixT, VectorT>(params);

    return Utopia::Finalize();
}

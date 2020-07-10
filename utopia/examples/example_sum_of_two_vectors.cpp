#include <iostream>
#include "utopia.hpp"

// Utopia may use different backend libraries
// to define a vector, a template class allows
// flexibility.
template <class Vector>
void sum_of_two_vectors() {
    using namespace utopia;

    // Initialisation of the communicator. The communicator
    // allows different processes to communicate between them,
    // sharing local information for their computation.
    using Comm = typename Traits<Vector>::Communicator;

    // SyzeType is a type that refers to n_locals and
    // n_global, which we will see in a second.
    using SizeType = typename Traits<Vector>::SizeType;

    // Scalar is a number which we can use to manipulate
    // a vector.
    using Scalar = typename Traits<Vector>::Scalar;

    // We declare the vectors.
    Vector a, b;

    // Default initialisation of the communicator.
    Comm& comm = Comm::get_default();

    // For each process we declare the number of local
    // entries that we want to use to compute our calculations.
    // In this case, our vector will be 10x1.
    SizeType n_local = 10;

    // We compute the global size of our vectors,
    // which we define as the number of process
    // multiplied by the number of local entries.
    SizeType n_global = n_local * comm.size();

    // We need a layout, that is a description of
    // the context in which the proccesses will
    // operate. A layout is therefore based on
    // the local entries involved in the task, the
    // communicator between them, and the global size of
    // the operation.
    auto l = layout(comm, n_local, n_global);

    // Initialise vector a. Note: this type of initialisation
    // assign the same value to each entry.
    // a will therefore be [5 5 5 5 5 5 5 5 5 5]'.
    a.values(l, 5);

    // Initialise vector b.
    b.values(l, 2);

    // ----- Initialising vectors with different values ------


    // We have to way to initialise vectors in a parallel way.
    // We can either do it locally or globally. A local initia-
    // -lisation involves a fixed range for the entries. Each
    // entry may have a different processor. On the other hand,
    // with a global intialisation we have the sum of all the
    // n_locals. 
    // You can run the programm to see the difference in terms
    // of output. The first output is local, the second is global.

    // In both case, we need to access the device: the device
    // is where we are saving our entries, such as 'host' or
    // the gpu. 
    // Locally, this can be done with, for example: 
    // auto a_view = local_view_device(a);.
    // Globally, with auto a_view = view_device(a);.

    // In both cases, we use the set(index, value) function 
    // to set a value at a certain index. 
    // We can also use get(index) to get the value at certain 
    // index. 
    // the disp function allows to visualise the content of a
    // vector or a matrix. 

    // Use 'make -j4 complete' to compile. 
    // Then run with, for example, 4 processor with
    // /mpirun  -n 4 ./examples/your_example'


    {
        auto a_view = local_view_device(a);
        parallel_for(
            local_range_device(a), UTOPIA_LAMBDA(const SizeType& i) {
                const auto ai = a_view.get(i);
                a_view.set(i, i);
            });

    }

    // Display the content. 
    auto c_expr = a + b;
    Vector c = c_expr;

    disp(c);

    {
        auto a_view = view_device(a);
        parallel_for(
            range_device(a), UTOPIA_LAMBDA(const SizeType& i) {
                const Scalar val = i;
                // const Scalar val = (i == 0) ? 1e-14 : ((i < n / 2.0) ? -i : i);
                a_view.set(i, i);
            });
    }

    // We return the sum of the two vectors.
    auto d_expr = a + b;
    Vector d = d_expr;

    disp(d);


    // ---------------------------------------------------------------------
    // -------------- Wrong wasy to initialise a vector --------------------


    // This initialisation does not work either since it is also not performed 
    // in parallel. 
    // a_view.set(0, 1);
    // a_view.set(1, 10);
    // a_view.set(2, 33);


    // This does not work either since it is also not performed in 
    // parallel. 
    auto j = 0;
    for (auto i = 0; i < n_local; ++i) {
        a_view.set(i, j);
        j++;
    }



    // ---------------------------------------------------------------------


}

int main(int argc, char** argv) {
    using namespace utopia;

// Checking which backend is active.
#ifdef WITH_PETSC
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

    Utopia::Init(argc, argv);

    sum_of_two_vectors<VectorT>();

    return Utopia::Finalize();
}


#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

// Create fowrad declerations file..
// Seperate files(per class)...
// Keep names of functions.
//#include "utopia_Layout.hpp"

//#include "../core/interfaces/utopia_AbstractMatrix.hpp"

namespace utopia {
    // Forward declarations
    template <typename Scalar, typename SizeType>
    class AbstractVector;

    template <typename Scalar, typename SizeType>
    class AbstractMatrix;

    template <typename Scalar, typename SizeType>
    class AlgebraFactory;

    template <class Communicator, int Order, typename LocalSizeType, typename SizeType>
    class Layout;

    template <class SelfCommunicator, int Order, typename LocalSizeType, typename SizeType>
    class Layout;

    class Communicator;

    class SelfCommunicator;

#ifdef WITH_PETSC
    UTOPIA_FACTORY_REGISTER_VECTOR(PetscVector);
    UTOPIA_FACTORY_REGISTER_MATRIX(PetscMatrix);
#endif  // WITH_PETSC

#ifdef WITH_TRILINOS
    UTOPIA_FACTORY_REGISTER_VECTOR(TpetraVector);
#endif  // WITH_PETSC

}  // namespace utopia

namespace scripting {

    using Scalar = double;
    using SizeType = int;
    using Factory = utopia::AlgebraFactory<Scalar, SizeType>;
    using Communicator = class Communicator;
    using SelfCommunicator = class SelfCommunicator;
    using LocalSizeType = int;

    void init();
    void finalize();
    void print_info();
    // void create_serial_layout();

    class SparseMatrix {
    public:
        using MatrixImpl = utopia::AbstractMatrix<Scalar, SizeType>;

        SparseMatrix();
        ~SparseMatrix();
        void print_info();
        // void create_square_matrix(const SizeType& /*size*/);
        void create_identity_matrix(const SizeType&);

        // Utility functions -----
        bool is_empty();
        void clear_matrix();

        // Display functions----
        void disp();

        // Question...: Create special functions that create special matrices?

    private:
        MatrixImpl* impl_;
    };

    class Vector {
    public:
        using VectorImpl = utopia::AbstractVector<Scalar, SizeType>;

        Vector();
        ~Vector();
        void print_info();
        void create_vector(const SizeType& /*size*/, const Scalar& /*value*/);

        // Mutators-------
        // void set_value(const SizeType&, const Scalar&);

        // Utility functions -----
        bool is_empty();
        void clear_vector();

        // Display function-------
        void disp();

        // Vector sum(const Vector&);

    private:
        VectorImpl* impl_;
    };

    class Layout {
    public:
        using LayoutImpl = utopia::Layout<Communicator, 1, LocalSizeType, SizeType>;
        using LayoutImplSelf = utopia::Layout<SelfCommunicator, 1, LocalSizeType, SizeType>;

        // TODO:: figure out how to implement comm for parallel layout
        // void create_parallel_layout(Comm, const SizeType&, const SizeType&);

    private:
        LayoutImpl* impl_;
        LayoutImplSelf* implSelf_;
    };

    class Communicator {
    public:
        using CommImpl = utopia::Communicator;
        Communicator();
        ~Communicator();

    private:
        CommImpl* impl_;
    };

    class SelfCommunicator {
    public:
        using SelfCommImpl = utopia::SelfCommunicator;
        SelfCommunicator();
        ~SelfCommunicator();

    private:
        SelfCommImpl* impl_;
    };

}  // namespace scripting

#endif  // UTOPIA_SCRIPT_HPP

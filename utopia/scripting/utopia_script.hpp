#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

// #include "utopia_AbstractVector.hpp"

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
    using LocalSizeType = int;
    using Communicator = class Communicator;

    void init();
    void finalize();
    void print_info();

    class SparseMatrix {
    public:
        using MatrixImpl = utopia::AbstractMatrix<Scalar, SizeType>;

        SparseMatrix();
        ~SparseMatrix();
        void print_info();
        void clear();
        void disp();
        bool empty();

    private:
        MatrixImpl* impl_;
    };

    class Vector {
    public:
        using VectorImpl = utopia::AbstractVector<Scalar, SizeType>;

        Vector();
        ~Vector();
        void print_info();
        void clear();
        bool empty();
        void disp();
        void create_vector(const SizeType& /*size*/, const Scalar& /*value*/);

    private:
        VectorImpl* impl_;
    };

    class Layout {
    public:
        using LayoutImpl = utopia::Layout<Communicator, 1, LocalSizeType, SizeType>;

        Layout();
        ~Layout();

    private:
        LayoutImpl* impl_;
    };

    class Communicator {
    public:
        using CommunicatorImpl = utopia::Communicator;

        Communicator();
        ~Communicator();

    private:
        CommunicatorImpl* impl_;
    };

    class SelfCommunicator {
    public:
        using SelfCommunicatorImpl_ = utopia::Layout<Communicator, 1, LocalSizeType, SizeType>;

        SelfCommunicator();
        ~SelfCommunicator();

    private:
        SelfCommunicator* impl_;
    };

}  // namespace scripting

#endif  // UTOPIA_SCRIPT_HPP

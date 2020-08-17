#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

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
    using SelfCommunicator = class SelfCommunicator;
    using Layout = class Layout;

    void init();
    void finalize();
    void print_info();

    class SparseMatrix {
    public:
        using MatrixImpl = utopia::AbstractMatrix<Scalar, SizeType>;

        SparseMatrix();
        ~SparseMatrix();

        void disp();
        void set(const SizeType&, const SizeType&, const Scalar&);

        void print_info();
        void clear();
        bool empty();

    private:
        MatrixImpl* impl_;
    };

    class Vector {
    public:
        using VectorImpl = utopia::AbstractVector<Scalar, SizeType>;

        Vector();
        ~Vector();

        void disp();
        void set(const SizeType& /*position*/, const Scalar& /*value*/);
        Scalar get(const SizeType& /*position*/);
        void describe();
        bool empty();
        void clear();
        void create_vector(const SizeType& /*size*/, const Scalar& /*value*/);
        Scalar norm1();
        // void values(const Layout&, const Scalar&);

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
        int rank();
        int size();

    private:
        CommunicatorImpl* impl_;
    };

    class SelfCommunicator {
    public:
        using SelfCommunicatorImpl_ = utopia::SelfCommunicator;

        SelfCommunicator();
        ~SelfCommunicator();

        void get_default();
        int rank();
        int size();

    private:
        SelfCommunicator* impl_;
    };

}  // namespace scripting

#endif  // UTOPIA_SCRIPT_HPP

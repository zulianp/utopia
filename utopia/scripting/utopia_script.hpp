// #include <stdlib.h>
// #include <vector>

#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

//#include "utopia_AbstractVector.hpp"

void print_array(double *seq, int n);

namespace utopia {
    // Forward declarations

    template <typename Scalar, typename SizeType>
    class AbstractVector;

    template <typename Scalar, typename SizeType>
    class AbstractMatrix;

    template <typename Scalar, typename SizeType>
    class AlgebraFactory;

    template <class Comm, int Order, typename LocalSizeType_, typename SizeType_>
    class Layout;

    class Communicator;

    template <typename Derived, int Order>
    class Tensor;

    template <class T, int Order>
    class DeviceView;

}  // namespace utopia

namespace scripting {
    using Scalar = double;
    using SizeType = int;
    using LocalSizeType = int;
    using Factory = utopia::AlgebraFactory<Scalar, SizeType>;
    using Derived = utopia::AbstractVector<Scalar, SizeType>;
    //  using T = utopia::AbstractVector<Scalar, SizeType>;

    void init();
    void finalize();
    void print_info();

    class Communicator {
    public:
        using CommunicatorImpl = utopia::Communicator;

        Communicator();
        ~Communicator();

        CommunicatorImpl *get_communicator() const { return impl_; }

    private:
        CommunicatorImpl *impl_;
    };

    class Layout {
    public:
        using LayoutImpl = const utopia::Layout<utopia::Communicator, 1, LocalSizeType, SizeType>;

        Layout(const Communicator &comm, LocalSizeType local_size, SizeType global_size);
        ~Layout();

        LayoutImpl *get_layout() const { return impl_; }

    private:
        const LayoutImpl *impl_;
    };

    class Vector {
    public:
        using VectorImpl = utopia::AbstractVector<Scalar, SizeType>;

        Vector();
        ~Vector();
        void print_info();
        SizeType local_size();
        void values(const Layout &l, const Scalar &value);
        void add(const SizeType &i, const Scalar &value);
        void axpy(Scalar alpha, Vector *x);
        Scalar dot(const Vector *x) const;
        void describe() const;
        bool equals(const Vector *other, const Scalar tol) const;
        void set(const SizeType &i, const Scalar &value);
        void print_array(double *seq, int n);
        // void convert_into_uvector(std::vector<double> values, const Layout &l);
        void serial_uconversion(double *seq, int n);
        void parallel_uconversion(float *values, const Layout &l);
        void scale(const Scalar &a);
        // void double *from_utopia_to_carray();
        void write_into_carray(double *double_array);
        Scalar get(const SizeType &i) const;
        // void convert(int *seq, int n);

    private:
        VectorImpl *impl_;
    };

    class SparseMatrix {
    public:
        using MatrixImpl = utopia::AbstractMatrix<Scalar, SizeType>;

        SparseMatrix();
        ~SparseMatrix();
        void print_info();

    private:
        MatrixImpl *impl_;
    };

}  // namespace scripting

#endif  // UTOPIA_SCRIPT_HPP

#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

// #include "utopia_AbstractVector.hpp"

void print_array(double* seq, int n);

namespace utopia {
    // Forward declarations
    template <typename Scalar, typename SizeType>
    class AbstractVector;

    template <typename Scalar, typename SizeType>
    class AbstractMatrix;

    template <typename Scalar, typename SizeType>
    class AlgebraFactory;
}  // namespace utopia

namespace scripting {

    using Scalar = double;
    using SizeType = int;
    using Factory = utopia::AlgebraFactory<Scalar, SizeType>;

    void init();
    void finalize();
    void print_info();

    class SparseMatrix {
    public:
        using MatrixImpl = utopia::AbstractMatrix<Scalar, SizeType>;

        SparseMatrix();
        ~SparseMatrix();
        void print_info();

    private:
        MatrixImpl* impl_;
    };

    class Vector {
    public:
        using VectorImpl = utopia::AbstractVector<Scalar, SizeType>;

        Vector();
        ~Vector();
        void print_info();

    private:
        VectorImpl* impl_;
    };

}  // namespace scripting

#endif  // UTOPIA_SCRIPT_HPP

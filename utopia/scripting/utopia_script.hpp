#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

// #include "utopia_AbstractVector.hpp"

namespace utopia {
    //Forward declarations
    template<typename Scalar, typename SizeType>
    class AbstractVector;

    template<typename Scalar, typename SizeType>
    class AbstractMatrix;

    template<typename Scalar, typename SizeType>
    class AlgebraFactory;
}

namespace scripting {
    using Scalar   = double;
    using SizeType = int;

    void init();
    void finalize();
    void print_info();

    class SparseMatrix {
    public:

        using MatrixImpl = utopia::AbstractMatrix<Scalar, SizeType>;
        using Factory    = utopia::AlgebraFactory<Scalar, SizeType>;

        SparseMatrix();
        ~SparseMatrix();
        void print_info();
    private:
        MatrixImpl * impl_;
    };

    class Vector {
    public:
        using VectorImpl = utopia::AbstractVector<Scalar, SizeType>;
        using Factory    = utopia::AlgebraFactory<Scalar, SizeType>;

        Vector();
        ~Vector();
        void print_info();
    private:
        VectorImpl * impl_;
    };

}

#endif //UTOPIA_SCRIPT_HPP

#ifndef UTOPIA_SCRIPT_HPP
#define UTOPIA_SCRIPT_HPP

namespace utopia {
    //Forward declarations
    template<typename Scalar, typename SizeType>
    class AbstractVector;

    template<typename Scalar, typename SizeType>
    class AbstractMatrix;

    template<typename Scalar, typename SizeType>
    class AlgebraFactory;
}

// #include <memory>

namespace scripting {
    // void init(int argc, char *argv[]);

    void init();
    void finalize();
    void print_info();

    class SparseMatrix {
    public:
        using Scalar   = double;
        using SizeType = int;
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
        using Scalar   = double;
        using SizeType = int;
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

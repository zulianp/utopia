#ifndef UTOPIA_UTOPIA_BACKENDSPEC_HPP
#define UTOPIA_UTOPIA_BACKENDSPEC_HPP

namespace utopia {

    template<class Matrix, class Vector, class Scalar, class SizeType>
    class BackendSpecBase {
    public:

        //Your <Backend> class has to implement this
//    static Backend &Instance() {
//        static Backend instance;
//        return instance;
//    }

        bool assign(Vector & /*left*/, Vector &/*right*/) {
            assert(false); //TODO
            return false;
        }

        bool assign(Vector &/*left*/, const Vector &/*right*/) {
            assert(false); //TODO
            return false;
        }

        bool assign(Matrix &/*left*/, const Matrix &/*right*/) {
            assert(false); //TODO
            return false;
        }

        bool assign(Matrix &left, Matrix &right) {
            using std::move;
            left = move(right);
            return false;
        }

        inline static Range range(const Vector &/*v*/) {
            assert(false); //TODO
            return Range(0);
        }

        inline static Range rowRange(const Matrix &/*m*/) {
            assert(false); //TODO
            return Range(0);
        }

        inline static Range colRange(const Matrix &/*m*/) {
            assert(false); //TODO
            return Range(0);
        }

        bool assignTransposed(Matrix &/*left*/, const Matrix &/*right*/) {
            assert(false); //TODO
            return false;
        }

        bool assignFromRange(Matrix &/*left*/, const Matrix &/*right*/, const Range &/*rowRange*/, const Range &/*colRange*/) {
            assert(false); //TODO
            return false;
        }

        bool assignFromRange(Vector &/*left*/, const Vector &/*right*/, const Range &/*rowRange*/, const Range &/*colRange*/) {
            assert(false); //TODO
            return false;
        }

        bool assignToRange(Matrix &/*left*/, const Matrix &/*right*/, const Range &/*rowRange*/, const Range &/*colRange*/) {
            assert(false); //TODO
            return false;
        }

        bool assignToRange(Matrix &/*left*/, const Identity &/**/, const Range &/*rowRange*/, const Range &/*colRange*/) {
            assert(false); //TODO
            return false;
        }

        bool assignToRange(Vector &/*left*/, const Vector &/*right*/, const Range &/*rowRange*/, const Range &/*colRange*/) {
            assert(false); //TODO
            return false;
        }

        inline bool set(Vector &/*vec*/, const SizeType /* index*/, const Scalar /* value*/) {
            assert(false); //TODO
            return false;
        }

        inline Scalar get(const Vector &/*vec*/, const SizeType /*index*/) {
            assert(false); //TODO
            return 0;
        }

        inline Scalar get(const Matrix &/*mat*/, const SizeType /*row*/, const SizeType /*col*/) {
            assert(false); //TODO
            return 0;
        }

        inline bool set(Matrix & /*mat*/, const SizeType /* row*/, const SizeType /* col*/, const Scalar /* value*/) {
            assert(false); //TODO
            return false;
        }


        bool size(const Vector &/* v*/, Size & /* size*/) {
            assert(false); //TODO
            return false;
        }

        bool size(const Matrix &/* m*/, Size &/* size*/) {
            assert(false); //TODO
            return false;
        }

        bool build(Matrix &/* m*/, const Size &/* size*/, const Identity &) {
            assert(false); //TODO
            return false;
        }

        bool build(Matrix &/*m*/, const Size &/*size*/, const Values <Scalar> &/*values*/) {
            assert(false); //TODO
            return false;
        }

        bool build(Vector &/*v*/, const Size &/*size*/, const Values <Scalar> &/*values*/) {
            assert(false); //TODO
            return false;
        }


        //more generic overloads
        template<typename Ordinal>
        inline bool set(Vector &/*v*/, const std::vector <Ordinal> &/*indices*/, const std::vector <Scalar> &/*values*/) {
            assert(false); //TODO
            return false;
        }

        template<typename Ordinal>
        inline bool set(Matrix &/*m*/, const std::vector <Ordinal> &/*rows*/, const std::vector <Ordinal> &/*columns*/,
                        const std::vector <Scalar> &/*values*/) {

            assert(false); //TODO
            return false;
        }

        template<class Comparator>
        bool compare(const Vector & /*left*/, const Vector &/*right*/, const Comparator &/*comp*/) {
            assert(false); //TODO
            return false;
        }

        template<class Comparator>
        bool compare(const Matrix &/*left */, const Matrix &/*right*/, const Comparator &/*comp*/) {
            assert(false); //TODO
            return false;
        }
    };


    ///Vector-Vector Operations (also Matrix-Matrix where matrices are just considered as vectors)
    template<class Matrix, class Vector, class Scalar, class SizeType>
    class BackendSpecLevel1 {
    public:
        //////////////////////////// TENSOR ORDER 1 => VECTOR  ////////////////////////////
        bool apply(const Vector &/*left */, const Vector &/*right*/, const Minus &, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }

        bool apply(const Vector &/*left*/, const Vector &/*right*/, const Plus &, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }

        /// Vector and Scalar
        bool apply(const Vector &/*vec*/, const Scalar /*value*/, const Multiplies &, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }


        bool scal(const Scalar /*scale_factor*/, const Vector &/*v*/, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }

        //result = scale_factor * vec + result
        bool axpy(const Scalar /*scale_factor*/, const Vector &/*vec*/, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }

        //result = scale_factor*left + right
        bool zaxpy(const Scalar /* scale_factor*/, const Vector &/*left*/, const Vector &/*right*/, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }

        bool apply(const Scalar /*value*/, const Vector &/*vec*/, const Multiplies &, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }

        // Reductions
        Scalar reduce(const Vector &/*vec*/, const Plus &) {
            assert(false); //TODO
            return 0;
        }

        ///.... other associative operators e.g., Minus, Multiplies, ....

        Scalar dot(const Vector &/*left*/, const Vector &/*right*/) {
            assert(false); //TODO
            return 0;
        }

        Scalar norm2(const Vector /*vector*/) {
            assert(false); //TODO
            return 0;
        }

        //////////////////////////// TENSOR ORDER 2 => MATRIX  ////////////////////////////
        bool apply(const Matrix &/*left*/, const Matrix &/*right*/, const Plus &, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        bool apply(const Matrix &/*left*/, const Matrix &/*right*/, const Multiplies &, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }

        bool apply(const Matrix &/*left*/, const Matrix &/*right*/, const Minus &, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        bool apply(const Matrix &/*mat*/, const Scalar &/*value*/, const Multiplies &, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        bool apply(const Scalar &/*value*/, const Matrix &/*mat*/, const Multiplies &, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        /// Matrix and Scalar
        bool scal(const Scalar /*scale_factor*/, const Matrix &/*m*/, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        //result = scale_factor * vec + result
        bool axpy(const Scalar /*scale_factor*/, const Matrix &/*vec*/, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        //result = scale_factor*left + right
        bool zaxpy(const Scalar /*scale_factor*/, const Matrix &/*left*/, const Matrix &/*right*/, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        /// You can implement everything with the TENSOR Versions (generic) if the wrapped library permits
        /// Generic versions of the operations can be implemented as follows (for example)
        template<class Tensor, class Operation>
        bool apply(const Tensor &/*left*/, const Tensor &/*right*/, const Operation &, Tensor &/*result*/) {
            assert(false); //TODO
            return false;
        }

        template<class Tensor, class Operation>
        bool apply(const Tensor &/*left*/, const Tensor &/*right*/, Tensor &/*result*/, const Operation &/*op*/) {
            assert(false); //TODO
            return false;
        }

        //result = scale_factor * vec + result
        template<class Left, class Right>
        bool axpy(const Scalar /*scale_factor*/, const Left &/*vec*/, Right &/*result*/) {
            assert(false); //TODO
            return false;
        }


        //result = scale_factor*left + right
        template<class Left, class Right, class Result>
        bool zaxpy(const Scalar /*scale_factor*/, const Left &/*left*/, const Right &/*right*/, Result &/*result*/) {
            assert(false); //TODO
            return false;
        }
    };

    //FIXME use the generalized version when possible
    ///Matrix-Vector operations
    template<class Matrix, class Vector, class Scalar, class SizeType>
    class BackendSpecLevel2 {
    public:
        bool gemv(const Scalar /*alpha*/, const Matrix &/*left*/, const Vector &/*right*/, const Scalar /* beta*/, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }

        //mv
        bool apply(const Matrix &/*left*/, const Vector &/*right*/, const Multiplies &, Vector &/*result*/) {
            assert(false); //TODO
            return false;
        }
    };

    //Matrix-Matrix operations
    template<class Matrix, class Vector, class Scalar, class SizeType>
    class BackendSpecLevel3 {
    public:
        bool gemm(const Scalar /*alpha*/, const Matrix &/*left*/, const Matrix &/*right*/, const Scalar /*beta*/, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        bool gemm(const Scalar /* alpha*/, const Matrix & /* left*/, const Matrix & /* right*/, bool/* transpose_left*/,
                  bool /*transpose_right*/, const Scalar/* beta*/, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }

        bool gemm(const Scalar  /*alpha*/, const Vector &/*left*/, const Matrix &/*right*/, bool /*transpose_left*/,
                  bool /*transpose_right*/, const Scalar/* beta*/, Matrix &/*result*/) {
            assert(false); //TODO
            return false;
        }
    };
}
#endif //UTOPIA_UTOPIA_BACKENDSPEC_HPP

#ifndef UTOPIA_DYNAMIC_TYPE_DISTRIBUTED_MATRIX_HPP
#define UTOPIA_DYNAMIC_TYPE_DISTRIBUTED_MATRIX_HPP

#include "utopia_ForwardDeclarations.hpp"

namespace utopia {

    /**
     * TODO MOVE TO OWN FILE
     * @brief utility class for using dynamically typed tensor instead of concrete one
     */
    template<
        typename Scalar,
        typename SizeType,
        class DerivedMatrixType,
        class ConcreteVectorType
    >
    class DynamicTypeDistributedMatrix :
        public DistributedMatrix<Scalar, SizeType>,
        //concrete interfaces (static polymorphism)
        public BLAS1Tensor<DerivedMatrixType>,
        public BLAS2Matrix<DerivedMatrixType, ConcreteVectorType>,
        public BLAS3Matrix<DerivedMatrixType>
        {
    public:
        using DistributedMatrix = utopia::DistributedMatrix<Scalar, SizeType>;
        using DistributedVector = utopia::DistributedVector<Scalar, SizeType>;

        virtual ~DynamicTypeDistributedMatrix() {}

        //route abstract type methods to concrete methods
        void swap(DistributedMatrix &x) override
        {
            derived_().swap(static_cast<DerivedMatrixType &>(x));
        }

        void copy(const DistributedMatrix &x) override
        {
            derived_().copy(static_cast<const DerivedMatrixType &>(x));
        }

        void axpy(const Scalar &a, const DistributedMatrix &x) override
        {
            derived_().axpy(a, static_cast<const DerivedMatrixType &>(x));
        }

        // virtual void gemv(const bool transpose, const Scalar &alpha, const Vector &x, const Scalar
        void transpose(DistributedMatrix &C) const override
        {
            derived_().transpose(static_cast<DerivedMatrixType &>(C));
        }

        void multiply(const Scalar &alpha, const DistributedMatrix &B, DistributedMatrix &C) const override
        {
            derived_().axpy(
                alpha,
                static_cast<const DerivedMatrixType &>(B),
                static_cast<const DerivedMatrixType &>(C)
            );
        }

        Scalar dot(const DistributedMatrix &x) const override
        {
            return derived_().dot(static_cast<const DerivedMatrixType &>(x));
        }

        void gemv(const bool transpose, const Scalar &alpha, const DistributedVector &x, const Scalar &beta, DistributedVector &y) const override
        {
            derived_().gemv(
                transpose,
                alpha,
                static_cast<const ConcreteVectorType &>(x),
                beta,
                static_cast<ConcreteVectorType &>(y)
            );
        }

        /// C := alpha * op(A) * op(B)
        void multiply(
            const bool transpose_A,
            const bool transpose_B,
            const DistributedMatrix &B,
            DistributedMatrix &C) const override
        {
            derived_().multiply(
                transpose_A,
                transpose_B,
                static_cast<const DerivedMatrixType &>(B),
                static_cast<DerivedMatrixType &>(C)
            );
        }

    private:
        DerivedMatrixType &derived_()
        {
            return static_cast<DerivedMatrixType &>(*this);
        }

        const DerivedMatrixType &derived_() const
        {
            return static_cast<const DerivedMatrixType &>(*this);
        }

    };

}

#endif //UTOPIA_DYNAMIC_TYPE_DISTRIBUTED_MATRIX_HPP

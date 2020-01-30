#ifndef UTOPIA_ABSTRACT_MATRIX_HPP
#define UTOPIA_ABSTRACT_MATRIX_HPP

#include <iostream>

#include "utopia_Matrix.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Communicator.hpp"
#include "utopia_Operator.hpp"
#include "utopia_AbstractVector.hpp"

namespace utopia {

    template<typename Scalar_, typename SizeType_>
    class Traits< AbstractMatrix<Scalar_, SizeType_> > {
    public:
        using Scalar = Scalar_;
        using SizeType = SizeType_;
    };

    //parallel types, collective operations
    template<typename Scalar_, typename SizeType_>
    class AbstractMatrix :
    public DistributedMatrix<Scalar_, SizeType_>,
    // public PolymorphicMatrix,
    public SparseConstructible<Scalar_, SizeType_>,
    public Normed<Scalar_>,
    public Reducible<Scalar_>,
    public ReducibleMatrix<Scalar_, SizeType_>,
    public Transformable<Scalar_>,
    // // Static polymorphic types
    public BLAS1Tensor<AbstractMatrix<Scalar_, SizeType_>>,
    public BLAS2Matrix<AbstractMatrix<Scalar_, SizeType_>, AbstractVector<Scalar_, SizeType_>>,
    public BLAS3Matrix<AbstractMatrix<Scalar_, SizeType_>>,
    public Comparable<AbstractMatrix<Scalar_, SizeType_>>,
    public Operator<AbstractVector<Scalar_, SizeType_>>
    // public Tensor<AbstractMatrix<Scalar_, SizeType_>, 2>,
    // public Selectable<AbstractMatrix<Scalar_, SizeType_>, 2>
    {
    public:
        virtual ~AbstractMatrix() {}
    };

    template<class Matrix>
    class Wrapper<Matrix, 2> : public AbstractMatrix<
                                        typename Traits<Matrix>::Scalar,
                                        typename Traits<Matrix>::SizeType
                                        >
    {
    public:
        using Scalar    = typename Traits<Matrix>::Scalar;
        using SizeType  = typename Traits<Matrix>::SizeType;
        using Vector    = typename Traits<Matrix>::Vector;
        using VectorWrapper  = utopia::Wrapper<Vector, 1>;

        using AbstractMatrix = utopia::AbstractMatrix<Scalar, SizeType>;
        using AbstractVector = utopia::AbstractVector<Scalar, SizeType>;

        template<class... Args>
        Wrapper(Args &&...args)
        : impl_(std::make_shared<Matrix>(std::forward<Args>(args)...))
        {}

        template<class... Args>
        void construct(Args &&...args)
        {
            impl_ = std::make_shared<Matrix>(std::forward<Args>(args)...);
        }

        //locks
        inline void read_lock() override
        {
            impl_->read_lock();
        }

        inline void write_lock(WriteMode mode) override
        {
            impl_->write_lock(mode);
        }

        inline void read_unlock() override
        {
            impl_->read_unlock();
        }

        inline void write_unlock(WriteMode mode) override
        {
            impl_->write_unlock(mode);
        }

        //basic mutators
        inline void set(const SizeType &i, const SizeType &j, const Scalar &value) override
        {
            impl_->set(i, j, value);
        }

        inline void add(const SizeType &i, const SizeType &j, const Scalar &value) override
        {
            impl_->add(i, j, value);
        }

        //print function
        inline void describe() const override
        {
            impl_->describe();
        }

        //utility functions
        inline bool empty() const override
        {
            return impl_->empty();
        }

        inline void clear() override
        {
            return impl_->clear();
        }

        inline SizeType rows() const override
        {
            return impl_->rows();
        }

        inline SizeType cols() const override
        {
            return impl_->cols();
        }

        //basic collective mutators allowing to write on other processes (e.g. for FE assembly)
        inline void c_set(const SizeType &i, const SizeType &j, const Scalar &value) override
        {
            impl_->c_set(i, j, value);
        }

        inline void c_add(const SizeType &i, const SizeType &j, const Scalar &value) override
        {
            impl_->c_add(i, j, value);
        }

        inline Range row_range() const override
        {
            return impl_->row_range();
        }
        inline Range col_range() const override
        {
            return impl_->col_range();
        }

        inline SizeType local_rows() const override
        {
            return impl_->local_rows();
        }
        inline SizeType local_cols() const override
        {
            return impl_->local_cols();
        }

        inline Size local_size() const override
        {
            return impl_->local_size();
        }

        inline Size size() const override
        {
            return impl_->size();
        }

        inline Communicator &comm() override
        {
            return impl_->comm();
        }

        inline const Communicator &comm() const override
        {
            return impl_->comm();
        }

        /////////////////////////////////////////////

        inline Scalar norm_infty() const override
        {
            return impl_->norm_infty();
        }

        inline Scalar norm1() const override
        {
            return impl_->norm1();
        }

        inline Scalar norm2() const override
        {
            return impl_->norm2();
        }

        /////////////////////////////////////////////

        inline Scalar reduce(const Plus &op) const override
        {
            return impl_->reduce(op);
        }

        inline Scalar reduce(const Min &op)  const override
        {
            return impl_->reduce(op);
        }

        inline Scalar reduce(const Max &op)  const override
        {
            return impl_->reduce(op);
        }

        /////////////////////////////////////////////

        inline SizeType nnz(const Scalar tol = 0.0) const override
        {
            return impl_->nnz(tol);
        }

        inline Scalar trace() const override
        {
            return impl_->trace();
        }

        /////////////////////////////////////////////

        inline void transform(const Sqrt &op) override
        {
            return impl_->transform(op);
        }

        inline void transform(const Pow2 &op) override
        {
            return impl_->transform(op);
        }

        inline void transform(const Log &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Exp &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Cos &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Sin &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Abs &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Minus &op) override
        {
            return impl_->transform(op);
        }

        inline void transform(const Pow &op)  override
        {
            return impl_->transform(op);
        }

        inline void transform(const Reciprocal<Scalar> &op) override
        {
            return impl_->transform(op);
        }

        /////////////////////////////////////////////

        inline void swap(AbstractMatrix &x) override
        {
            auto &x_w = static_cast<const Wrapper &>(x);
            impl_->swap(*x_w.impl_);
        }

        ///<Scalar>SCAL - x = a*x
        inline void scale(const Scalar &a) override
        {
            impl_->scale(a);
        }

        ///<Scalar>COPY - copy x into y (this)
        inline void copy(const AbstractMatrix &x) override
        {
            auto &x_w = static_cast<const Wrapper &>(x);
            impl_->copy(*x_w.impl_);
        }

        ///<Scalar>AXPY - y = a*x + y
        inline void axpy(const Scalar &a, const AbstractMatrix &x) override
        {
            auto &x_w = static_cast<const Wrapper &>(x);
            impl_->axpy(a, *x_w.impl_);
        }

        ///<Scalar>DOT - dot product
        inline Scalar dot(const AbstractMatrix &x) const override
        {
            auto &x_w = static_cast<const Wrapper &>(x);
            return impl_->dot(*x_w.impl_);
        }

        /////////////////////////////////////////////


        /// y = A * x
        inline void multiply(const AbstractVector &x, AbstractVector &y) const override
        {
            auto &x_w = static_cast<const VectorWrapper &>(x);
            auto &y_w = static_cast<VectorWrapper &>(y);

            impl_->multiply(x_w.get(), y_w.get());
        }

        /// y = alpha * A * x
        inline void multiply(const Scalar &alpha, const AbstractVector &x, AbstractVector &y) const override
        {
            auto &x_w = static_cast<const VectorWrapper &>(x);
            auto &y_w = static_cast<VectorWrapper &>(y);

            impl_->multiply(alpha, x_w.get(), y_w.get());
        }

        /// y = alpha * A * x
        inline void transpose_multiply(const AbstractVector &x, AbstractVector &y) const override
        {
            auto &x_w = static_cast<const VectorWrapper &>(x);
            auto &y_w = static_cast<VectorWrapper &>(y);

            impl_->transpose_multiply(x_w.get(), y_w.get());
        }

        /// y = alpha * A^T * x
        inline void transpose_multiply(const Scalar &alpha, const AbstractVector &x, AbstractVector &y) const override
        {
            auto &x_w = static_cast<const VectorWrapper &>(x);
            auto &y_w = static_cast<VectorWrapper &>(y);

            impl_->transpose_multiply(alpha, x_w.get(), y_w.get());
        }

        /// y := alpha * A * x + beta * y
        inline void multiply_add(const Scalar &alpha, const AbstractVector &x, const Scalar &beta, AbstractVector &y) const override
        {
            auto &x_w = static_cast<const VectorWrapper &>(x);
            auto &y_w = static_cast<VectorWrapper &>(y);

            impl_->multiply_add(alpha, x_w.get(), beta, y_w.get());
        }

        /// y := alpha * A' * x + beta * y
        inline void transpose_multiply_add(const Scalar &alpha, const AbstractVector &x, const Scalar &beta, AbstractVector &y) const override
        {
            auto &x_w = static_cast<const VectorWrapper &>(x);
            auto &y_w = static_cast<VectorWrapper &>(y);

            impl_->transpose_multiply_add(alpha, x_w.get(), beta, y_w.get());
        }

        /// y := alpha * op(A) * x + beta * y
        inline void gemv(
            const bool transpose,
            const Scalar &alpha,
            const AbstractVector &x,
            const Scalar &beta,
            AbstractVector &y
        ) const override
        {
            auto &x_w = static_cast<const VectorWrapper &>(x);
            auto &y_w = static_cast<VectorWrapper &>(y);

            impl_->gemv(transpose, alpha, x_w.get(), beta, y_w.get());
        }

        /////////////////////////////////////////////


        inline void transpose(AbstractMatrix &C) const override
        {
            auto &C_w = static_cast<Wrapper &>(C);
            impl_->transpose(C_w.get());
        }

        inline void multiply(const AbstractMatrix &B, AbstractMatrix &C) const override
        {
            auto &B_w = static_cast<const Wrapper &>(B);
            auto &C_w = static_cast<Wrapper &>(C);
            impl_->multiply(B_w.get(), C_w.get());
        }

        /// C := alpha * A * B
        inline void multiply(const Scalar &alpha, const AbstractMatrix &B, AbstractMatrix &C) const override
        {
            auto &B_w = static_cast<const Wrapper &>(B);
            auto &C_w = static_cast<Wrapper &>(C);
            impl_->multiply(alpha, B_w.get(), C_w.get());
        }

        /// C := A^T * B
        inline void transpose_multiply(const AbstractMatrix &B, AbstractMatrix &C) const override
        {
            auto &B_w = static_cast<const Wrapper &>(B);
            auto &C_w = static_cast<Wrapper &>(C);
            impl_->transpose_multiply(B_w.get(), C_w.get());
        }

        /// C := A * B^T
        inline void multiply_transpose(const AbstractMatrix &B, AbstractMatrix &C) const override
        {
            auto &B_w = static_cast<const Wrapper &>(B);
            auto &C_w = static_cast<Wrapper &>(C);
            impl_->multiply_transpose(B_w.get(), C_w.get());
        }

        /// C := alpha * op(A) * op(B)
        inline void multiply(
            const bool transpose_A,
            const bool transpose_B,
            const AbstractMatrix &B,
            AbstractMatrix &C) const override
        {
            auto &B_w = static_cast<const Wrapper &>(B);
            auto &C_w = static_cast<Wrapper &>(C);
            impl_->multiply(transpose_A, transpose_B, B_w.get(), C_w.get());
        }

        /////////////////////////////////////////////

        inline bool equals(const AbstractMatrix &other, const Scalar &tol = 0.0) const override
        {
            auto &other_w = static_cast<const Wrapper &>(other);
            return impl_->equals(*other_w.impl_, tol);
        }

        /////////////////////////////////////////////

        inline bool apply(const AbstractVector &x, AbstractVector &y) const
        {
            auto &x_w = static_cast<const VectorWrapper &>(x);
            auto &y_w = static_cast<VectorWrapper &>(y);
            return impl_->apply(x_w.get(), y_w.get());
        }

        /////////////////////////////////////////////

        inline void sparse(const Size &s, const SizeType &nnz) override
        {
            impl_->sparse(s, nnz);
        }

        inline void local_sparse(const Size &s, const SizeType &nnz) override
        {
            impl_->local_sparse(s, nnz);
        }

        inline void identity(const Size &s, const Scalar &diag = 1.0) override
        {
            impl_->identity(s, diag);
        }

        inline void local_identity(const Size &s, const Scalar &diag = 1.0) override
        {
            impl_->local_identity(s, diag);
        }

        /////////////////////////////////////////////

        Matrix &get()
        {
            return *impl_;
        }

        const Matrix &get() const
        {
            return *impl_;
        }

        std::shared_ptr<Matrix> ptr()
        {
            return impl_;
        }

        std::shared_ptr<const Matrix> ptr() const
        {
            return impl_;
        }

    private:
        std::shared_ptr<Matrix> impl_;

    };
}

#endif //UTOPIA_ABSTRACT_MATRIX_HPP



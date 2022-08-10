#ifndef UTOPIA_UTOPIA_PETSCMATRIX_H
#define UTOPIA_UTOPIA_PETSCMATRIX_H

// Core includes
#include "utopia_Allocations.hpp"
#include "utopia_ArrayView.hpp"
#include "utopia_BLAS_Operands.hpp"
#include "utopia_Base.hpp"
#include "utopia_Comparable.hpp"
#include "utopia_Constructible.hpp"
#include "utopia_ElementWiseOperand.hpp"
#include "utopia_Matrix.hpp"
#include "utopia_Normed.hpp"
#include "utopia_Operator.hpp"
#include "utopia_Range.hpp"
#include "utopia_Reducible.hpp"
#include "utopia_Select.hpp"
#include "utopia_Size.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Tracer.hpp"
#include "utopia_Transformable.hpp"

// Backend includes
#include "utopia_DynamicTypeDistributedMatrix.hpp"
#include "utopia_petsc_Base.hpp"
#include "utopia_petsc_Communicator.hpp"
#include "utopia_petsc_Error.hpp"
#include "utopia_petsc_ForwardDeclarations.hpp"
#include "utopia_petsc_IndexSet.hpp"
#include "utopia_petsc_Traits.hpp"

#include <memory>

#include "petscmat.h"
#include "petscsys.h"

namespace utopia {
    class PetscMatrixMemory {
    public:
        PetscMatrixMemory(const MPI_Comm comm = PETSC_COMM_WORLD) { MatCreate(comm, &_mat); }

        PetscMatrixMemory(Mat &mat, const bool owner = false) : _mat(mat), owner_(owner) {}

        ~PetscMatrixMemory() { destroy(); }

        inline void destroy() {
            if (owner_) {
                MatDestroy(&_mat);
            }
        }

        inline Mat &implementation() { return _mat; }

        inline const Mat &implementation() const {
            assert(_mat != nullptr);
            return _mat;
        }

        // MAT_DO_NOT_COPY_VALUES or MAT_COPY_VALUES, cause it to copy the numerical values in the matrix
        // MAT_SHARE_NONZERO_PATTERN
        inline void duplicate(PetscMatrixMemory &other, MatDuplicateOption opt = MAT_COPY_VALUES) const {
            // if(other.owner_) {
            // MatDestroy(&other._mat);
            // }

            other.destroy();

            PetscErrorHandler::Check(MatDuplicate(_mat, opt, &other._mat));

            UTOPIA_REPORT_ALLOC("PetscMatrix::duplicate");
        }

        inline void convert(PetscMatrixMemory &other, MatType newtype) {
            // MAT_REUSE_MATRIX is only supported for inplace conversion, otherwise use MAT_INITIAL_MATRIX.
            PetscErrorHandler::Check(MatConvert(_mat, newtype, MAT_INITIAL_MATRIX, &other._mat));

            UTOPIA_REPORT_ALLOC("PetscMatrix::convert");
        }

        inline void convert(MatType newtype) {
            // MAT_REUSE_MATRIX is only supported for inplace conversion, otherwise use MAT_INITIAL_MATRIX.
            PetscErrorHandler::Check(MatConvert(_mat, newtype, MAT_REUSE_MATRIX, &_mat));
        }

        PetscMatrixMemory(const PetscMatrixMemory &other) {
            PetscErrorHandler::Check(MatCopy(other._mat, _mat, SAME_NONZERO_PATTERN));
        }

        inline MPI_Comm communicator() const {
            MPI_Comm comm = PetscObjectComm((PetscObject)implementation());
            assert(comm != MPI_COMM_NULL);
            return comm;
        }

        void set_owner(const bool owner) { owner_ = owner; }

        inline bool is_owner() const { return owner_; }

    private:
        Mat _mat{nullptr};
        bool owner_{true};
    };

    class PetscMatrix :
        // Dynamic polymorphic types
        public DistributedMatrix<PetscScalar, PetscInt>,
        public PolymorphicMatrix,
        public Normed<PetscScalar>,
        public Reducible<PetscScalar>,
        public ReducibleMatrix<PetscScalar, PetscInt>,
        public Transformable<PetscScalar>,
        // Static polymorphic types
        public Constructible<PetscMatrix>,
        public BLAS1Tensor<PetscMatrix>,
        public BLAS2Matrix<PetscMatrix, PetscVector>,
        public BLAS3Matrix<PetscMatrix>,
        public Comparable<PetscMatrix>,
        public Operator<PetscVector>,
        public Tensor<PetscMatrix, 2>,
        public Selectable<PetscMatrix, 2> {
    public:
        using Scalar = PetscScalar;
        using SizeType = PetscInt;
        using Super = utopia::Tensor<PetscMatrix, 2>;
        using Super::Super;
        using MatrixLayout = typename Traits<PetscMatrix>::MatrixLayout;
        using IndexArray = typename Traits<PetscMatrix>::IndexArray;
        using ScalarArray = typename Traits<PetscMatrix>::ScalarArray;

        ////////////////////////////////////////////////////////////////////
        ///////////////////////// BOILERPLATE CODE FOR EDSL ////////////////
        ////////////////////////////////////////////////////////////////////
        void init_empty(const PetscCommunicator &comm) { wrapper_ = std::make_shared<PetscMatrixMemory>(comm.get()); }

        // nvcc has problems with a default argument for PetscMatrix(const PetscCommunicator &comm), so we trick it by
        // forwarding the default constructor to the constructor with the default argument
        PetscMatrix() : PetscMatrix(PetscCommunicator::get_default()) {}
        PetscMatrix(const PetscCommunicator &comm /*= PetscCommunicator::get_default() */) : comm_(comm) {
            init_empty(comm);
        }

        /// same as sparse
        PetscMatrix(const MatrixLayout &layout, const SizeType &d_nnz, const SizeType &o_nnz) {
            init_empty(layout.comm());
            sparse(layout, d_nnz, o_nnz);
        }

        /// same as dense(l, 0.0)
        PetscMatrix(const MatrixLayout &layout) {
            init_empty(layout.comm());
            dense(layout, 0.0);
        }

        PetscMatrix(PetscMatrix &&other) : wrapper_(std::move(other.wrapper_)) { this->update_mirror(); }

        PetscMatrix(const PetscMatrix &other) {
            using std::make_shared;

            if (other.empty()) {
                wrapper_ = make_shared<PetscMatrixMemory>(other.communicator());
                return;
            }

            wrapper_ = make_shared<PetscMatrixMemory>();
            other.wrapper_->duplicate(*wrapper_);
            this->update_mirror();
        }

        PetscMatrix &operator=(const PetscMatrix &other) {
            if (wrapper_ == other.wrapper_) return *this;

            if (other.empty()) {
                clear();
                return *this;
            }

            wrapper_ = std::make_shared<PetscMatrixMemory>();
            other.wrapper_->duplicate(*wrapper_);
            this->update_mirror();
            return *this;
        }

        PetscMatrix &operator=(PetscMatrix &&other) {
            if (wrapper_ == other.wrapper_) return *this;

            this->wrapper_ = other.wrapper_;
            other.wrapper_ = nullptr;
            this->update_mirror();
            return *this;
        }

        template <class Expr>
        PetscMatrix(const Expression<Expr> &expr) {
            static_assert(!std::is_same<Expr, Tensor<PetscVector, 1>>::value, "cannot assign a vector to a matrix");

            // FIXME see if expression can provide a comm
            init_empty(comm_);
            // THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
            // assert(valid());
        }

        template <class Expr>
        inline PetscMatrix &operator=(const Expression<Expr> &expr) {
            Super::assign_eval(expr.derived());
            // assert(valid());
            return *this;
        }

        void assign(const PetscMatrix &other) override { copy(other); }

        void assign(PetscMatrix &&other) override {
            comm_ = std::move(other.comm_);
            wrapper_ = std::move(other.wrapper_);
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR DistributedObject ////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        PetscCommunicator &comm() override { return comm_; }

        const PetscCommunicator &comm() const override { return comm_; }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR MatrixBase and DistributedMatrix /////////////
        ///////////////////////////////////////////////////////////////////////////

        inline void c_set(const SizeType &i, const SizeType &j, const Scalar &value) override {
            check_error(MatSetValue(implementation(), i, j, value, INSERT_VALUES));
        }

        inline void c_set_block(const SizeType &i, const SizeType &j, const Scalar *block) {
            check_error(MatSetValuesBlocked(implementation(), 1, &i, 1, &j, block, INSERT_VALUES));
        }

        inline void c_add_block(const SizeType &i, const SizeType &j, const Scalar *block) {
            check_error(MatSetValuesBlocked(implementation(), 1, &i, 1, &j, block, ADD_VALUES));
        }

        inline void c_add(const SizeType &i, const SizeType &j, const Scalar &value) override {
            check_error(MatSetValue(implementation(), i, j, value, ADD_VALUES));
        }

        inline Range row_range() const override {
            SizeType r_begin, r_end;
            MatGetOwnershipRange(implementation(), &r_begin, &r_end);
            assert(Range(r_begin, r_end).valid());
            return {r_begin, r_end};
        }

        inline Range col_range() const override {
            SizeType r_begin, r_end;
            MatGetOwnershipRangeColumn(implementation(), &r_begin, &r_end);
            assert(Range(r_begin, r_end).valid());
            return {r_begin, r_end};
        }

        inline ArrayView<const PetscInt> col_ranges() const {
            const PetscInt *ranges;
            MatGetOwnershipRangesColumn(this->raw_type(), &ranges);
            return ArrayView<const PetscInt>(ranges, comm().size());
        }

        inline ArrayView<const PetscInt> row_ranges() const {
            const PetscInt *ranges;
            MatGetOwnershipRanges(this->raw_type(), &ranges);
            return ArrayView<const PetscInt>(ranges, comm().size());
        }

        inline SizeType rows() const override {
            SizeType ret;
            check_error(MatGetSize(implementation(), &ret, nullptr));
            return ret;
        }

        inline SizeType cols() const override {
            SizeType ret;
            check_error(MatGetSize(implementation(), nullptr, &ret));
            return ret;
        }

        inline SizeType local_rows() const override {
            SizeType ret;
            check_error(MatGetLocalSize(implementation(), &ret, nullptr));
            return ret;
        }

        inline SizeType local_cols() const override {
            SizeType ret;
            check_error(MatGetLocalSize(implementation(), nullptr, &ret));
            return ret;
        }

        inline Size size() const override {
            SizeType rows, cols;
            check_error(MatGetSize(implementation(), &rows, &cols));
            return {rows, cols};
        }

        inline Size local_size() const override {
            SizeType rows, cols;
            check_error(MatGetLocalSize(implementation(), &rows, &cols));
            return {rows, cols};
        }

        void clear() override;
        bool empty() const override;

        inline void set(const SizeType &row, const SizeType &col, const Scalar &value) override {
            assert(row_range().inside(row));
            check_error(MatSetValue(implementation(), row, col, value, INSERT_VALUES));
        }

        inline void add(const SizeType &row, const SizeType &col, const Scalar &value) override {
            assert(row_range().inside(row));
            check_error(MatSetValue(implementation(), row, col, value, ADD_VALUES));
        }

        inline void describe() const override {
            if (empty()) {
                utopia::out() << "PetscMatrix: Empty" << std::endl;
            }

            MatView(wrapper_->implementation(), PETSC_VIEWER_STDOUT_(communicator()));
        }

        // REMOVEME
        inline MPI_Comm communicator() const { return wrapper_->communicator(); }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Constructible ////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline void sparse(const MatrixLayout &layout,
                           const SizeType &nnz_d_block,
                           const SizeType &nnz_o_block) override {
            comm_ = layout.comm();
            matij_init(comm().get(),
                       MATAIJ,
                       layout.local_size(0),
                       layout.local_size(1),
                       layout.size(0),
                       layout.size(1),
                       nnz_d_block,
                       nnz_o_block);
        }

        inline void sparse(const MatrixLayout &layout, const IndexArray &d_nnz, const IndexArray &o_nnz) {
            comm_ = layout.comm();
            matij_init(comm().get(),
                       MATAIJ,
                       layout.local_size(0),
                       layout.local_size(1),
                       layout.size(0),
                       layout.size(1),
                       d_nnz,
                       o_nnz);
        }

        void crs(const MatrixLayout &layout,
                 const IndexArray &row_ptr,
                 const IndexArray &col_idx,
                 const ScalarArray &values);

        inline void block_sparse(const MatrixLayout &layout,
                                 const IndexArray &d_nnz,
                                 const IndexArray &o_nnz,
                                 const SizeType block_size) {
            comm_ = layout.comm();
            mat_baij_init(comm().get(),
                          layout.local_size(0),
                          layout.local_size(1),
                          layout.size(0),
                          layout.size(1),
                          d_nnz,
                          o_nnz,
                          block_size);
        }

        void identity(const Scalar &diag = 1.0);

        inline void identity(const MatrixLayout &layout, const Scalar &diag = 1.0) override {
            comm_ = layout.comm();
            matij_init_identity(
                comm().get(), MATAIJ, layout.local_size(0), layout.local_size(1), layout.size(0), layout.size(1), diag);
        }

        inline void dense(const MatrixLayout &layout, const Scalar &val = 0.0) override {
            comm_ = layout.comm();
            dense_init_values(comm().get(),
                              MATDENSE,
                              layout.local_size(0),
                              layout.local_size(1),
                              layout.size(0),
                              layout.size(1),
                              val);
        }

        inline void dense_identity(const MatrixLayout &layout, const Scalar &diag = 1.0) override {
            comm_ = layout.comm();
            dense_init_identity(comm().get(),
                                MATDENSE,
                                layout.local_size(0),
                                layout.local_size(1),
                                layout.size(0),
                                layout.size(1),
                                diag);
        }

        void convert_to_scalar_matrix();
        void convert_to_scalar_matrix(PetscMatrix &scalar_matrix);

        // void identity(const Scalar &diag = 1.0);

        // inline void identity(const Size &s, const Scalar &diag = 1.0) override {
        //     matij_init_identity(comm().get(), MATAIJ, PETSC_DECIDE, PETSC_DECIDE, s.get(0), s.get(1), diag);
        // }

        // inline void dense_identity(const Size &s, const Scalar &diag = 1.0) override {
        //     dense_init_identity(comm().get(), MATDENSE, PETSC_DECIDE, PETSC_DECIDE, s.get(0), s.get(1), diag);
        // }

        // inline void local_dense_identity(const Size &s, const Scalar &diag = 1.0) override {
        //     dense_init_identity(comm().get(), MATDENSE, s.get(0), s.get(1), PETSC_DETERMINE, PETSC_DETERMINE, diag);
        // }

        // inline void values(const Size &s, const Scalar &val) override {
        //     dense_init_values(comm().get(), MATDENSE, PETSC_DECIDE, PETSC_DECIDE, s.get(0), s.get(1), val);
        // }

        // inline void sparse(const Size &s, const SizeType &nnz) override {
        //     matij_init(comm().get(), MATAIJ, PETSC_DECIDE, PETSC_DECIDE, s.get(0), s.get(1), nnz, nnz);
        // }

        // inline void local_identity(const Size &s, const Scalar &diag = 1.0) override {
        //     matij_init_identity(comm().get(), MATAIJ, s.get(0), s.get(1), PETSC_DETERMINE, PETSC_DETERMINE, diag);
        // }

        // inline void local_values(const Size &s, const Scalar &val) override {
        //     dense_init_values(comm().get(), MATDENSE, s.get(0), s.get(1), PETSC_DETERMINE, PETSC_DETERMINE, val);
        // }

        // /// Specialize for sparse matrices
        // inline void local_sparse(const Size &s, const SizeType &nnz) override {
        //     matij_init(comm().get(), MATAIJ, s.get(0), s.get(1), PETSC_DETERMINE, PETSC_DETERMINE, nnz, nnz);
        // }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Normed ////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline Scalar norm1() const override {
            Scalar val;
            MatNorm(implementation(), NORM_1, &val);
            return val;
        }

        inline Scalar norm2() const override {
            Scalar val;
            MatNorm(implementation(), NORM_FROBENIUS, &val);
            return val;
        }

        inline Scalar norm_infty() const override {
            Scalar val;
            MatNorm(implementation(), NORM_INFINITY, &val);
            return val;
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Reducible ////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        Scalar sum() const override;
        Scalar max() const override;
        Scalar min() const override;

        inline Scalar reduce(const Plus &) const override { return sum(); }
        inline Scalar reduce(const Min &) const override { return min(); }
        inline Scalar reduce(const Max &) const override { return max(); }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR ReducbileMatrix //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        SizeType nnz(const Scalar tol = 0.0) const override;

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Transformable //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        void transform(const Sqrt &) override;
        void transform(const Pow2 &) override;
        void transform(const Log &) override;
        void transform(const Exp &) override;
        void transform(const Cos &) override;
        void transform(const Sin &) override;
        void transform(const Abs &) override;
        void transform(const Minus &) override;

        void transform(const Pow &p) override;
        void transform(const Reciprocal<Scalar> &op) override;

        void transform(std::function<Scalar(const Scalar &)> op);
        void transform(std::function<Scalar(const SizeType &, const SizeType &, const Scalar &)> op);

        // helper
        template <class F>
        void transform_values(F f);

        template <class F>
        void transform_ijv(F f);

        template <class Op>
        void read(Op op) const;

        template <class Op>
        void read_reverse(Op op) const;

        template <class Op, class MPIOp>
        Scalar parallel_reduce_values(const Op &op, const MPIOp &mpi_op, const Scalar &initial_value) const;

        template <class Map, class Reduce, class MPIOp, typename Accumulator>
        void map_reduce(const Map &map, const Reduce &reduce, const MPIOp &, Accumulator &accumulator) const;

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR ElementWiseOperand //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        // virtual void e_mul(const PetscMatrix &other) = 0;
        // virtual void e_min(const PetscMatrix &other) = 0;
        // virtual void e_max(const PetscMatrix &other) = 0;

        // FIXME inherit from ElementWiseOperand
        /// x = x*(1./a)
        inline void e_div(const Scalar &a)  // override
        {
            check_error(MatScale(implementation(), 1. / a));
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR BLAS1Tensor //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        ///< Scalar>SWAP - swap x and y
        inline void swap(PetscMatrix &x) override {
            using std::swap;
            swap(comm_, x.comm_);
            swap(wrapper_, x.wrapper_);
        }

        ///< Scalar>SCAL - x = a*x
        inline void scale(const Scalar &a) override { check_error(MatScale(implementation(), a)); }

        ///< Scalar>COPY - copy x into y (this)
        void copy(const PetscMatrix &other) override;

        ///< Scalar>COPY - copy x into y (this)
        void same_nnz_pattern_copy(const PetscMatrix &other);

        ///< Scalar>AXPY - y = a*x + y
        void axpy(const Scalar &a, const PetscMatrix &x) override;

        ///< Scalar>DOT - dot product
        Scalar dot(const PetscMatrix & /*other*/) const override {
            assert(false && "IMPLEMENT ME");
            return -1.0;
        }

        /// I<Scalar>AMAX - index of max abs value
        // SizeType amax() const override
        // {
        //    // MatGetRowMaxAbs()
        //    assert(false && "IMPLEMENT ME");
        //    return 1.0;
        // }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Operator //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline bool apply(const PetscVector &in, PetscVector &out) const override {
            if (empty()) return false;

            this->multiply(in, out);

            return true;
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR BLAS2Matrix //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        //////////////////
        ///// IN SUPPORT OF OVERRIDES
        /// result = v2 + this * v1
        void multiply_add(const PetscVector &v1, const PetscVector &v2, PetscVector &result) const;
        void transpose_multiply_add(const PetscVector &v1, const PetscVector &v2, PetscVector &result) const;
        //////////////////

        /// y = A * x
        void multiply(const PetscVector &vec, PetscVector &result) const override;

        /// y = A * x
        void transpose_multiply(const PetscVector &vec, PetscVector &result) const override;

        /// y = alpha * A * x
        void multiply(const Scalar &alpha, const PetscVector &x, PetscVector &y) const override;

        /// y = alpha * A^T * x
        void transpose_multiply(const Scalar &alpha, const PetscVector &x, PetscVector &y) const override;

        /// y := alpha * A * x + beta * y
        void multiply_add(const Scalar &alpha, const PetscVector &x, const Scalar &beta, PetscVector &y) const override;

        /// y := alpha * A' * x + beta * y
        void transpose_multiply_add(const Scalar &alpha,
                                    const PetscVector &x,
                                    const Scalar &beta,
                                    PetscVector &y) const override;

        void gemv(const bool transpose,
                  const Scalar &alpha,
                  const PetscVector &x,
                  const Scalar &beta,
                  PetscVector &y) const override;

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR BLAS3Matrix //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        void transpose();
        void transpose(PetscMatrix &result) const override;

        void multiply(const PetscMatrix &mat, PetscMatrix &result) const override;

        /// C := alpha * A * B
        void multiply(const Scalar &alpha, const PetscMatrix &B, PetscMatrix &C) const override;

        /// C := alpha * A^T * B
        void transpose_multiply(const PetscMatrix &mat, PetscMatrix &result) const override;

        /// C := alpha * A * B^T
        void multiply_transpose(const PetscMatrix &mat, PetscMatrix &result) const override;

        /// C := alpha * op(A) * op(B)
        void multiply(const bool transpose_A,
                      const bool transpose_B,
                      const PetscMatrix &B,
                      PetscMatrix &C) const override;

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Comparable //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        bool equals(const PetscMatrix &other, const Scalar &tol = 0.0) const override;

        ////////////////////////////////////////////////////////////////////

        Mat &implementation() { return wrapper_->implementation(); }

        inline MatType type() const {
            MatType type;
            MatGetType(implementation(), &type);
            return type;
        }

        virtual MatType type_override() const;

        inline std::string name() const {
            const char *name;
            PetscObjectGetName((PetscObject)implementation(), &name);
            return name;
        }

        inline void rename(const std::string &name) { PetscObjectSetName((PetscObject)implementation(), name.c_str()); }

        void wrap(Mat &mat) { wrapper_ = std::make_shared<PetscMatrixMemory>(mat, false); }

        const Mat &implementation() const { return wrapper_->implementation(); }

        Mat &raw_type() { return wrapper_->implementation(); }

        const Mat &raw_type() const { return wrapper_->implementation(); }

        inline Scalar get(const SizeType row, const SizeType col) const {
            Scalar value;
            check_error(MatGetValues(implementation(), 1, &row, 1, &col, &value));
            return value;
        }

        template <class Index, class Values>
        void add_matrix(const Index &rows, const Index &cols, const Values &values) {
            // assert(rows.size() * cols.size() == values.size());

            check_error(MatSetValues(raw_type(),
                                     static_cast<PetscInt>(rows.size()),
                                     &rows[0],
                                     static_cast<PetscInt>(cols.size()),
                                     &cols[0],
                                     &values[0],
                                     ADD_VALUES));
        }

        template <class Index, class Values>
        void add_matrix_blocked(const Index &rows, const Index &cols, const Values &values) {
            // assert(rows.size() * cols.size() == values.size());

            check_error(MatSetValuesBlocked(raw_type(),
                                            static_cast<PetscInt>(rows.size()),
                                            &rows[0],
                                            static_cast<PetscInt>(cols.size()),
                                            &cols[0],
                                            &values[0],
                                            ADD_VALUES));
        }

        void add_matrix(const std::vector<SizeType> &rows,
                        const std::vector<SizeType> &cols,
                        const std::vector<Scalar> &values);

        void set_matrix(const std::vector<SizeType> &rows,
                        const std::vector<SizeType> &cols,
                        const std::vector<Scalar> &values);

        inline void write_lock(WriteMode) override {}

        inline void write_unlock(WriteMode) override {
            UTOPIA_TRACE_REGION_BEGIN("PetscMatrix::write_unlock(...)");

            check_error(MatAssemblyBegin(implementation(), MAT_FINAL_ASSEMBLY));
            check_error(MatAssemblyEnd(implementation(), MAT_FINAL_ASSEMBLY));

            UTOPIA_TRACE_REGION_END("PetscMatrix::write_unlock(...)");
        }

        inline void read_lock() override {}
        inline void read_unlock() override {}

        void dense_init(MPI_Comm comm,
                        MatType dense_type,
                        SizeType rows_local,
                        SizeType cols_local,
                        SizeType rows_global,
                        SizeType cols_global);

        bool write(const std::string &path) const;
        bool write_binary(const std::string &path) const;
        bool write_matlab(const std::string &path) const;

        inline bool read(const std::string &path) { return read(comm().get(), path); }

        bool read(MPI_Comm comm, const std::string &path);
        void copy_from(Mat mat);
        void copy_to(Mat mat) const;
        void copy_to(Mat *mat) const;

        bool is_sparse() const override;

        void select(const PetscIndexSet &row_index, const PetscIndexSet &col_index, PetscMatrix &result) const override;

        void local_select(const Range &local_row_range,
                          const Range &local_col_range,
                          const Range &global_col_range,
                          PetscMatrix &result) const;

        void select(const Range &global_row_range, const Range &global_col_range, PetscMatrix &result) const;

        inline Scalar trace() const override {
            Scalar ret;
            check_error(MatGetTrace(implementation(), &ret));
            return ret;
        }

        void row_sum(PetscVector &col) const;
        void row_max(PetscVector &col) const;
        void row_abs_max(PetscVector &col) const;
        void row_min(PetscVector &col) const;

        void col_sum(PetscVector &col) const;
        void col_max(PetscVector &) const { assert(false && "IMPLEMENT ME"); }
        void col_min(PetscVector &) const { assert(false && "IMPLEMENT ME"); }

        bool is_mpi() const;
        bool has_nan_or_inf() const;

        inline bool is_owner() const {
            if (wrapper_) {
                return wrapper_->is_owner();
            } else {
                return false;
            }
        }

        // petsc says that it is correct only for square matrices
        void build_diag(PetscVector &result) const;
        void build_diag(PetscMatrix &result) const;
        void diag(const PetscVector &other);
        void diag(const PetscMatrix &other) { other.build_diag(*this); }

        void col(const SizeType id, PetscVector &result) const;

        inline void shift_diag(const Scalar factor) { check_error(MatShift(implementation(), factor)); }
        void set_diag(const PetscVector &d);

        void shift_diag(const PetscVector &d);
        void dense_init_diag(MatType dense_type, const PetscVector &diag);
        void matij_init_diag(const PetscVector &diag);

        void dense_init_values(MPI_Comm comm,
                               MatType dense_type,
                               SizeType local_rows,
                               SizeType local_cols,
                               SizeType global_rows,
                               SizeType global_cols,
                               Scalar value);

        void dense_init_identity(MPI_Comm comm,
                                 MatType dense_type,
                                 SizeType local_rows,
                                 SizeType local_cols,
                                 SizeType global_rows,
                                 SizeType global_cols,
                                 Scalar scale_factor);

        void matij_init_identity(MPI_Comm comm,
                                 MatType sparse_type,
                                 SizeType local_rows,
                                 SizeType local_cols,
                                 SizeType global_rows,
                                 SizeType global_cols,
                                 Scalar scale_factor);

        void matij_init(MPI_Comm comm,
                        MatType type,
                        SizeType rows_local,
                        SizeType cols_local,
                        SizeType rows_global,
                        SizeType cols_global,
                        SizeType d_nnz,
                        SizeType o_nnz);

        void matij_init(MPI_Comm comm,
                        MatType type,
                        SizeType rows_local,
                        SizeType cols_local,
                        SizeType rows_global,
                        SizeType cols_global,
                        const IndexArray &d_nnz,
                        const IndexArray &o_nnz);

        template <class Int>
        void matij_init(MPI_Comm comm,
                        MatType type,
                        SizeType rows_local,
                        SizeType cols_local,
                        SizeType rows_global,
                        SizeType cols_global,
                        const std::vector<Int> &d_nnz,
                        const std::vector<Int> &o_nnz) {
            std::vector<SizeType> petsc_d_nnz(d_nnz.size()), petsc_o_nnz(o_nnz.size());

            std::copy(d_nnz.begin(), d_nnz.end(), petsc_d_nnz.begin());
            std::copy(o_nnz.begin(), o_nnz.end(), petsc_o_nnz.begin());

            matij_init(comm, type, rows_local, cols_local, rows_global, cols_global, petsc_d_nnz, petsc_o_nnz);
        }

        void mat_baij_init(MPI_Comm comm,
                           SizeType rows_local,
                           SizeType cols_local,
                           SizeType rows_global,
                           SizeType cols_global,
                           SizeType d_nnz,
                           SizeType o_nnz,
                           SizeType block_size);

        void mat_baij_init(MPI_Comm comm,
                           SizeType rows_local,
                           SizeType cols_local,
                           SizeType rows_global,
                           SizeType cols_global,
                           const IndexArray &d_nnz,
                           const IndexArray &o_nnz,
                           SizeType block_size);

        void nest(MPI_Comm comm,
                  SizeType nr,
                  const IS is_row[],
                  SizeType nc,
                  const IS is_col[],
                  const Mat a[],
                  const bool use_mat_nest_type = false);

        // void mat_aij_cusparse_init(
        //  	MPI_Comm comm,
        //  	SizeType rows_local,
        //  	SizeType cols_local,
        //  	SizeType rows_global,
        //  	SizeType cols_global,
        //  	SizeType d_nnz,
        //  	SizeType o_nnz
        // );

        inline void destroy() { wrapper_->destroy(); }

        void inverse(PetscMatrix &result) const;

        void convert_to_mat_baij(const PetscInt block_size, PetscMatrix &output);
        void convert_to_mat_baij(const SizeType block_size);

        bool is_initialized_as(MPI_Comm comm,
                               MatType dense_type,
                               SizeType local_rows,
                               SizeType local_cols,
                               SizeType global_rows,
                               SizeType global_cols);

        bool has_type(MatType type) const;
        bool same_type(const PetscMatrix &other) const;
        bool is_cuda() const;
        bool is_assembled() const;

        static bool is_block(Mat mat);
        bool is_block() const;

        void convert_from(const Mat &mat);
        void convert_to(Mat &mat) const;

        void set_zero_rows(const PetscIndexSet &idx, const Scalar &diag = 0.0);

        // only for dense
        void set(const Scalar &value);

        void diagonal_block(PetscMatrix &other) const;

        void diag_scale_right(const PetscVector &diag);
        void diag_scale_left(const PetscVector &diag);

        inline std::string get_class() const override { return "PetscMatrix"; }

        void assign(const Range &row_range, const Range &col_range, const PetscMatrix &block);

        SizeType global_nnz() const;
        SizeType local_nnz() const;

        inline bool is_alias(const PetscMatrix &other) const {
            if (empty()) return false;
            if (other.empty()) return false;
            return raw_type() == other.raw_type();
        }

        void update_mirror();

    private:
        PetscCommunicator comm_;
        std::shared_ptr<PetscMatrixMemory> wrapper_;

        inline static bool check_error(const SizeType err) { return PetscErrorHandler::Check(err); }

        void select_aux(const std::vector<SizeType> &row_index,
                        const std::vector<SizeType> &col_index,
                        PetscMatrix &result) const;

        void par_assign_from_local_is(const std::vector<SizeType> &remote_rows,
                                      const std::vector<SizeType> &remote_cols,
                                      const SizeType global_col_offset,
                                      const Range &local_col_range,
                                      PetscMatrix &result) const;

        VecType compatible_cuda_vec_type() const;

        // y = A * x;
        bool create_vecs(Vec *x, Vec *y) const;

        bool valid() const;

        // helpers
        template <class Op>
        void op_transform(const Op &op);

        template <class F>
        void transform_values_seqaij(F f);

        template <class F>
        void transform_ijv_seqaij(F f);

        template <class F>
        void transform_values_mpiaij(F op);

        template <class F>
        void transform_ijv_mpiaij(F f);

        // wait for petsc version
        // template<class F>
        // void read_ijv_seqaij(F op);
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_PETSCMATRIX_H

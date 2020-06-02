#ifndef UTOPIA_TPETRAMATRIX_H
#define UTOPIA_TPETRAMATRIX_H

#include "utopia_Base.hpp"
#include "utopia_Logger.hpp"
#include "utopia_Matrix.hpp"
#include "utopia_Normed.hpp"
#include "utopia_Operator.hpp"
#include "utopia_Range.hpp"
#include "utopia_Select.hpp"
#include "utopia_Size.hpp"
#include "utopia_Tensor.hpp"
#include "utopia_Traits.hpp"

#include "utopia_Tpetra_Vector.hpp"
#include "utopia_trilinos_Communicator.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include "utopia_DynamicTypeDistributedMatrix.hpp"

#include <iostream>
#include <memory>

namespace utopia {
    // template<class NodeType>
    class TpetraMatrix
        : public DistributedSparseMatrix<TpetraScalar, TpetraSizeType>,
          public SparseConstructible<TpetraScalar, TpetraSizeType>,
          public Normed<TpetraScalar>,
          public Transformable<TpetraScalar>,
          // Static polymorphic types
          // public DynamicTypeDistributedMatrix<TpetraScalar, TpetraSizeType, TpetraMatrix, TpetraVector>,
          public BLAS1Tensor<TpetraMatrix>,
          public BLAS2Matrix<TpetraMatrix, TpetraVector>,
          public BLAS3Matrix<TpetraMatrix>,
          public Comparable<TpetraMatrix>,
          public Operator<TpetraVector>,
          public Tensor<TpetraMatrix, 2>,
          public Selectable<TpetraMatrix, 2> {
    public:
        /////////////////////////////////////////////////////////////
        // typedef definitions
        /////////////////////////////////////////////////////////////

        using Scalar = Traits<TpetraMatrix>::Scalar;
        using SizeType = Traits<TpetraMatrix>::SizeType;
        using LocalSizeType = Traits<TpetraMatrix>::LocalSizeType;
        using Node = Traits<TpetraMatrix>::Node;
        using IndexSet = Traits<TpetraMatrix>::IndexSet;
        using IndexArray = Traits<TpetraMatrix>::IndexArray;
        using ScalarArray = Traits<TpetraMatrix>::ScalarArray;
        using MatrixLayout = Traits<TpetraMatrix>::MatrixLayout;

        using BLAS2Matrix<TpetraMatrix, TpetraVector>::multiply;
        using BLAS2Matrix<TpetraMatrix, TpetraVector>::transpose_multiply;

        // types of Trilinos Objects
        using CrsMatrixType = Tpetra::CrsMatrix<Scalar, LocalSizeType, SizeType, Node>;
        using RCPCrsMatrixType = Teuchos::RCP<CrsMatrixType>;
        using RCPCommType = Teuchos::RCP<const Teuchos::Comm<int>>;
        using MapType = Tpetra::Map<LocalSizeType, SizeType, Node>;
        using RCPMapType = Teuchos::RCP<const MapType>;

        /////////////////////////////////////////////////////////////
        // Constructors
        /////////////////////////////////////////////////////////////

        // Default Constructor
        TpetraMatrix() : owner_(true) {}

        // deep copy
        //     template <class Node2>
        // RCPCrsMatrixType  clone (
        //   const Teuchos::RCP<Node2> & node2,
        //   const Teuchos::RCP<Teuchos::ParameterList> & params = Teuchos::null)
        // {
        //   return mat_->clone(node2, params);
        // }

        TpetraMatrix(const TpetraMatrix &other) : owner_(true) {
            if (!other.is_null()) {
                mat_.reset(new CrsMatrixType(other.implementation(), Teuchos::Copy));
            }
        }

        TpetraMatrix(TpetraMatrix &&other) : mat_(std::move(other.mat_)), owner_(std::move(other.owner_)) {}

        TpetraMatrix(const RCPCrsMatrixType &mat, const bool owner = false) : mat_(mat), owner_(owner) {}

        /////////////////////////////////////////////////////////////
        // Destructor
        /////////////////////////////////////////////////////////////

        ~TpetraMatrix() override {}

        ////////////////////////////////////////////////////////////////////
        ///////////////////////// BOILERPLATE CODE FOR EDSL ////////////////
        ////////////////////////////////////////////////////////////////////

        using Super = utopia::Tensor<TpetraMatrix, 2>;
        using Super::Super;

        template <class Expr>
        TpetraMatrix(const Expression<Expr> &expr) : owner_(true) {
            static_assert(!std::is_same<Expr, Tensor<TpetraVector, 1>>::value, "cannot assign a vector to a matrix");
            // THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
        }

        template <class Expr>
        inline TpetraMatrix &operator=(const Expression<Expr> &expr) {
            Super::assign_eval(expr.derived());
            return *this;
        }

        void assign(const TpetraMatrix &other) override { copy(other); }

        void assign(TpetraMatrix &&other) override {
            owner_ = std::move(other.owner_);
            comm_ = std::move(other.comm_);
            mat_ = std::move(other.mat_);
        }

        void copy(const TpetraMatrix &other) override {
            if (this == &other) return;

            if (other.is_null()) {
                mat_.reset();
                owner_ = true;
                return;
            }

            mat_.reset(new CrsMatrixType(other.implementation(), Teuchos::Copy));
            owner_ = true;
        }

        void select(const IndexSet &row_index, const IndexSet &col_index, TpetraMatrix &result) const override;

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR DistributedObject ////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        TrilinosCommunicator &comm() override { return comm_; }

        const TrilinosCommunicator &comm() const override { return comm_; }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR MatrixBase and DistributedMatrix /////////////
        ///////////////////////////////////////////////////////////////////////////

        void c_set(const SizeType &row, const SizeType &col, const Scalar &value) override;
        void c_add(const SizeType &row, const SizeType &col, const Scalar &value) override;

        SizeType rows() const override;
        SizeType cols() const override;
        SizeType local_rows() const override;
        SizeType local_cols() const override;

        inline Range row_range() const override {
            return {implementation().getRowMap()->getMinGlobalIndex(),
                    implementation().getRowMap()->getMaxGlobalIndex() + 1};
        }

        inline Range col_range() const override {
            if (implementation().getDomainMap().is_null()) {
                assert(!init_->domain_map.is_null());
                return {init_->domain_map->getMinGlobalIndex(), init_->domain_map->getMaxGlobalIndex() + 1};
            } else {
                return {implementation().getDomainMap()->getMinGlobalIndex(),
                        implementation().getDomainMap()->getMaxGlobalIndex() + 1};
            }
        }

        inline Size size() const override {
            if (is_null()) {
                return {0, 0};
            }

            if (implementation().isFillComplete()) {
                return {implementation().getGlobalNumRows(), implementation().getGlobalNumCols()};
            } else {
                assert(!implementation().getRowMap().is_null());

                if (implementation().getDomainMap().is_null()) {
                    assert(!init_->domain_map.is_null());
                    return {implementation().getRowMap()->getGlobalNumElements(),
                            init_->domain_map->getGlobalNumElements()};
                } else {
                    return {implementation().getRowMap()->getGlobalNumElements(),
                            implementation().getDomainMap()->getGlobalNumElements()};
                }
            }
        }

        inline bool empty() const override { return is_null(); }

        inline Size local_size() const override {
            if (is_null()) {
                return {0, 0};
            }

            assert(!implementation().getRowMap().is_null());

            if (implementation().getDomainMap().is_null()) {
                assert(!init_->domain_map.is_null());
                return {implementation().getRowMap()->getNodeNumElements(), init_->domain_map->getNodeNumElements()};
            } else {
                return {implementation().getRowMap()->getNodeNumElements(),
                        implementation().getDomainMap()->getNodeNumElements()};
            }
        }

        void describe(std::ostream &os) const {
            auto out = Teuchos::getFancyOStream(Teuchos::rcpFromRef(os));
            implementation().describe(*out, Teuchos::EVerbosityLevel::VERB_EXTREME);
        }

        inline void describe() const override { describe(std::cout); }

        void clear() override;

        /////////////////////////////////////////////////////////////
        // OVERRIDES for SparseConstructible
        /////////////////////////////////////////////////////////////

        inline void sparse(const MatrixLayout &layout, const SizeType nnz_d_block, const SizeType nnz_o_block) {
            comm_ = layout.comm();
            crs_init(comm_.get(),
                     layout.local_size(0),
                     layout.local_size(1),
                     layout.size(0),
                     layout.size(1),
                     std::max(nnz_d_block, nnz_o_block));
        }

        inline void identity(const MatrixLayout &layout, const Scalar &diag = 1.0) {
            comm_ = layout.comm();
            crs_identity(comm_.get(), layout.local_size(0), layout.local_size(1), layout.size(0), layout.size(1), diag);
        }

        void identity(const Size &s, const Scalar &diag = 1.0) override;

        /// Specialize for sparse matrices
        void sparse(const Size &s, const SizeType & /*nnz*/) override;

        /// Specialize for sparse matrices
        void local_sparse(const Size &s, const SizeType & /*nnz*/) override;

        void local_identity(const Size &s, const Scalar &diag = 1.0) override;

        /////////////////////////////////////////////////////////////
        // OVERRIDES for Normed
        /////////////////////////////////////////////////////////////

        Scalar norm_infty() const override;
        Scalar norm1() const override;
        Scalar norm2() const override;

        /////////////////////////////////////////////////////////////
        // Overloading Operators
        /////////////////////////////////////////////////////////////

        TpetraMatrix &operator=(const TpetraMatrix &other) {
            copy(other);
            return *this;
        }

        TpetraMatrix &operator=(TpetraMatrix &&other) {
            if (this == &other) return *this;

            if (other.is_null()) {
                mat_.reset();
                owner_ = true;
                return *this;
            }

            mat_ = std::move(other.mat_);
            owner_ = std::move(other.owner_);
            return *this;
        }

        void finalize();

        // void fillComplete()
        // {
        //     mat_->fillComplete();
        // }

        // void replaceGlobalValues (const SizeType globalRow, const LocalSizeType numEnt, const SC vals[], const
        // SizeType cols[])
        // {
        //     mat_->replaceGlobalValues(globalRow, numEnt, vals, cols);
        // }

        // void replaceLocalValues (const LocalSizeType localRow, const LocalSizeType numEnt,  const SC vals[], const
        // LocalSizeType cols[] )
        // {
        //     mat_->replaceLocalValues(localRow, numEnt, vals, cols);
        // }

        RCPCommType communicator() const { return implementation().getMap()->getComm(); }

        void set_owner(const bool owner) { owner_ = owner; }

        // API functions
        void crs_init(const RCPCommType &comm,
                      std::size_t rows_local,
                      std::size_t cols_local,
                      Tpetra::global_size_t rows_global,
                      Tpetra::global_size_t cols_global,
                      std::size_t nnz_x_row);

        void crs_init(const RCPCommType &comm,
                      std::size_t rows_local,
                      std::size_t cols_local,
                      Tpetra::global_size_t rows_global,
                      Tpetra::global_size_t cols_global,
                      const Teuchos::ArrayRCP<size_t> &rowPtr,
                      const Teuchos::ArrayRCP<LocalSizeType> &cols,
                      const Teuchos::ArrayRCP<Scalar> &values);

        void crs_identity(const RCPCommType &comm,
                          std::size_t rows_local,
                          std::size_t cols_local,
                          Tpetra::global_size_t rows_global,
                          Tpetra::global_size_t cols_global,
                          const Scalar factor);

        void crs(const MatrixLayout &l,
                 const Teuchos::ArrayRCP<size_t> &rowPtr,
                 const Teuchos::ArrayRCP<LocalSizeType> &cols,
                 const Teuchos::ArrayRCP<Scalar> &values) {
            comm_ = l.comm();

            crs_init(comm_.get(), l.local_size(0), l.local_size(1), l.size(0), l.size(1), rowPtr, cols, values);
        }

        inline void read_lock() override {
            // TODO?
        }

        inline void read_unlock() override {
            // TODO?
        }

        inline void write_lock(WriteMode mode = utopia::AUTO) override {
            UTOPIA_UNUSED(mode);
            // TODO?
            implementation().resumeFill();
        }

        inline void write_unlock(WriteMode mode = utopia::AUTO) override {
            UTOPIA_UNUSED(mode);
            this->finalize();
        }

        void set(const SizeType &row, const SizeType &col, const Scalar &value) override;
        Scalar get(const SizeType &row, const SizeType &col) const;
        void add(const SizeType &row, const SizeType &col, const Scalar &value) override;

        template <typename Integer>
        void add_matrix(const std::vector<Integer> &rows,
                        const std::vector<Integer> &cols,
                        const std::vector<Scalar> &values) {
            // FIXME and find more efficient way
            const auto n_rows = rows.size();
            const auto n_cols = cols.size();
            assert(values.size() == n_rows * n_cols);

            for (std::size_t i = 0; i < n_rows; ++i) {
                const auto i_offset = i * n_rows;
                for (std::size_t j = 0; j < n_cols; ++j) {
                    add(rows[i], cols[j], values[i_offset + j]);
                }
            }
        }

        template <typename Integer>
        void set_matrix(const std::vector<Integer> &rows,
                        const std::vector<Integer> &cols,
                        const std::vector<Scalar> &values) {
            // FIXME and find more efficient way
            const auto n_rows = rows.size();
            const auto n_cols = cols.size();
            assert(values.size() == n_rows * n_cols);

            for (std::size_t i = 0; i < n_rows; ++i) {
                const auto i_offset = i * n_rows;
                for (std::size_t j = 0; j < n_cols; ++j) {
                    set(rows[i], cols[j], values[i_offset + j]);
                }
            }
        }

        template <typename Integer>
        void set_matrix(const std::vector<Integer> &rows, const std::vector<Integer> &cols, const Scalar value) {
            // FIXME and find more efficient way
            const auto n_rows = rows.size();
            const auto n_cols = cols.size();
            for (std::size_t i = 0; i < n_rows; ++i) {
                for (std::size_t j = 0; j < n_cols; ++j) {
                    set(rows[i], cols[j], value);
                }
            }
        }

        inline bool apply(const TpetraVector &in, TpetraVector &out) const override {
            if (empty()) return false;

            this->multiply(in, out);

            return true;
        }

        void multiply(const TpetraVector &vec, TpetraVector &result) const override;
        void transpose_multiply(const TpetraVector &vec, TpetraVector &result) const override;

        void multiply(const TpetraMatrix &right, TpetraMatrix &result) const override;
        void multiply(const Scalar &alpha, const TpetraMatrix &B, TpetraMatrix &C) const override;
        // result = tranpose(*this) * mat
        void transpose_multiply(const TpetraMatrix &right, TpetraMatrix &result) const override;

        void multiply(const bool transpose_this,
                      const bool transpose_right,
                      const TpetraMatrix &right,
                      TpetraMatrix &result) const override;

        void axpy(const Scalar &alpha, const TpetraMatrix &x) override;
        void transpose(TpetraMatrix &mat) const override;

        void build_diag(TpetraVector &d) const;
        void diag(const TpetraVector &d);
        void diag(const TpetraMatrix &mat);

        void diag_scale_left(const TpetraVector &d);

        inline RCPCrsMatrixType &raw_type() { return implementation_ptr(); }

        inline const RCPCrsMatrixType &raw_type() const { return implementation_ptr(); }

        inline CrsMatrixType &implementation() {
            assert(!mat_.is_null());
            return *mat_;
        }

        inline const CrsMatrixType &implementation() const {
            assert(!mat_.is_null());
            return *mat_;
        }

        inline RCPCrsMatrixType &implementation_ptr() { return mat_; }

        inline const RCPCrsMatrixType &implementation_ptr() const {
            assert(!mat_.is_null());
            return mat_;
        }

        inline bool is_null() const { return mat_.is_null(); }

        inline bool read(const std::string &path) { return read(comm().get(), path); }

        bool read(const Teuchos::RCP<const Teuchos::Comm<int>> &comm, const std::string &path);
        bool write(const std::string &path) const;

        bool is_valid(const bool verbose = false) const;

        Scalar sum() const;
        // FIXME
        inline Scalar reduce(const Plus &) const  // override
        {
            return sum();
        }

        void set_domain_and_range(const RCPMapType &domain_map, const RCPMapType &range_map) {
            init_ = std::make_shared<InitStructs>();
            init_->domain_map = domain_map;
            init_->range_map = range_map;
        }

        void build_from_structure(const TpetraMatrix &rhs);

        void set_zero_rows(const IndexSet &index, const Scalar &diag);

        ////////////////////////////////////////////////////////////////////////
        //////////////////////////// Transormable //////////////////////////////
        ////////////////////////////////////////////////////////////////////////

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

        template <class Op>
        Scalar local_parallel_reduce_values(Op op, const Scalar &initial_value) const;

        template <class Op, class MPIOp>
        Scalar parallel_reduce_values(Op op, MPIOp mpi_op, const Scalar &initial_value) const;

        template <class F>
        void transform_values(F op);

        template <class Op>
        void transform_ijv(Op op);

        template <class F>
        void read(F op) const;

        ////////////////////////////////////////////////////////////////////////
        //////////////////////////// Blas1Tensor //////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        void swap(TpetraMatrix &x) override;
        void scale(const Scalar &alpha) override;
        Scalar dot(const TpetraMatrix &other) const override;

        ////////////////////////////////////////////////////////////////////////
        //////////////////////////// Blas2Matrix //////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        void gemv(const bool transpose,
                  const Scalar &alpha,
                  const TpetraVector &x,
                  const Scalar &beta,
                  TpetraVector &y) const override;

        ////////////////////////////////////////////////////////////////////////
        //////////////////////////// Comparable //////////////////////////////
        ////////////////////////////////////////////////////////////////////////

        bool equals(const TpetraMatrix &other, const Scalar &tol = 0.0) const override;

        inline std::string get_class() const override { return "TpetraMatrix"; }

        inline bool is_alias(const TpetraMatrix &other) const { return mat_ == other.mat_; }

        void shift_diag(const TpetraVector &d);

    private:
        TrilinosCommunicator comm_;
        RCPCrsMatrixType mat_;
        bool owner_;

        typedef struct {
            RCPMapType domain_map;
            RCPMapType range_map;
        } InitStructs;

        std::shared_ptr<InitStructs> init_;

    NVCC_PRIVATE
        template <class Op>
        void aux_transform(const Op &op);

    };  // TpetraMatrix
}  // namespace utopia

#endif  // UTOPIA_TPETRAMATRIX_H

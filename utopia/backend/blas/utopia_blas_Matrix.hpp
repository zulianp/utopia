#ifndef UTOPIA_BLAS_MATRIX_HPP
#define UTOPIA_BLAS_MATRIX_HPP

#include "utopia_Interfaces.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_blas_Traits.hpp"
#include "utopia_BLAS_Operands.hpp"
#include "utopia_blas_Algorithms.hpp"
#include "utopia_blas_ForwardDeclarations.hpp"
#include "utopia_blas_Vector.hpp"
#include "utopia_Size.hpp"
#include "utopia_Constructible.hpp"
#include "utopia_ElementWiseOperand.hpp"
#include "utopia_Normed.hpp"
#include "utopia_Transformable.hpp"
#include "utopia_Reducible.hpp"
#include "utopia_blas_IndexSet.hpp"


#include <vector>
#include <memory>
#include <iostream>

namespace utopia {
    template<typename T>
    class BlasDenseMatrix : 
        // Dynamic polymorphic types
        public DenseMatrix<T, std::size_t>,
        public ReducibleMatrix<T, std::size_t>,
        public Constructible<T, std::size_t, 2>,
        public Normed<T>,
        public Transformable<T>,
        public Reducible<T>,
        // Static polymorphic types
        public Tensor<BlasDenseMatrix<T>, 2>,
        public BLAS1Tensor<BlasDenseMatrix<T>>,
        public BLAS2Matrix<BlasDenseMatrix<T>, BlasVector<T>>,
        public BLAS3Matrix<BlasDenseMatrix<T>>,
        public Comparable<BlasDenseMatrix<T>>,
        public ElementWiseOperand<BlasDenseMatrix<T>> {
    public:
        typedef std::vector<T> Entries;
        using SizeType = std::size_t;
        using Scalar = T;
       
        using BlasVector = utopia::BlasVector<T>;
        using BLAS2Matrix<BlasDenseMatrix, BlasVector>::multiply;
        using BLAS3Matrix<BlasDenseMatrix>::multiply;
        using BLAS2Matrix<BlasDenseMatrix, BlasVector>::transpose_multiply;
        using BLAS3Matrix<BlasDenseMatrix>::transpose_multiply;


        ////////////////////////////////////////////////////////////////////
        ///////////////////////// BOILERPLATE CODE FOR EDSL ////////////////
        ////////////////////////////////////////////////////////////////////

        using Super = utopia::Tensor<BlasDenseMatrix<T>, 2>;

        template<class Expr>
        BlasDenseMatrix(const Expression<Expr> &expr)
        : rows_(0), cols_(0)
        {
            //THIS HAS TO BE HERE IN EVERY UTOPIA TENSOR CLASS
            Super::construct_eval(expr.derived());
        }

        template<class Expr>
        inline BlasDenseMatrix &operator=(const Expression<Expr> &expr)
        {
            Super::assign_eval(expr.derived());
            return *this;
        }

        void assign(const BlasDenseMatrix &other) override
        {
            copy(other);
        }

        void assign(BlasDenseMatrix &&other) override
        {
            entries_ = std::move(other.entries_);
            rows_ = std::move(other.rows_);
            cols_ = std::move(other.cols_);
        }

        ////////////////////////////////////////////////////////////////////

        BlasDenseMatrix(BlasDenseMatrix &&other)
        : entries_(std::move(other.entries_)), 
          rows_(std::move(other.rows_)),
          cols_(std::move(other.cols_))
        {}

        BlasDenseMatrix(const BlasDenseMatrix &other)
        : entries_(other.entries_), 
          rows_(other.rows_),
          cols_(other.cols_)
        {}

        inline BlasDenseMatrix &operator=(const BlasDenseMatrix &other)
        {
            if(this == &other) return *this;
            copy(other);
            return *this;
        }

        inline BlasDenseMatrix &operator=(BlasDenseMatrix &&other)
        {
            if(this == &other) return *this;
            entries_ = std::move(other.entries_);
            rows_    = std::move(other.rows_);
            cols_    = std::move(other.cols_);
            return *this;
        }

        BlasDenseMatrix(SizeType rows, SizeType cols) : entries_(rows * cols), rows_(rows), cols_(cols)
        {}

        BlasDenseMatrix(SizeType rows, SizeType cols, T value) : entries_(rows * cols, value), rows_(rows), cols_(cols)
        {}

        ~BlasDenseMatrix() { }

        BlasDenseMatrix() : rows_(0), cols_(0) {}

        BlasDenseMatrix(SizeType rows, SizeType cols, std::initializer_list<T> args) : rows_(rows), cols_(cols)
        {
            using std::copy;

            entries_.resize(rows_*cols_);
            copy(args.begin(), args.end(), entries_.begin());
        }

        BlasDenseMatrix(const Entries& e)
        {
            rows_ = e.n_elements();
            cols_ = 1;
            entries_ = e;
        }

        SizeType rows() const override {
            return rows_;
        }

        SizeType cols() const override {
            return cols_;
        }

        SizeType n_elements() const {
            return rows() * cols();
        }

        const Entries &entries() const {
            return entries_;
        }

        Entries &entries() {
            return entries_;
        }

        template<class EntriesT>
        void set_entries(EntriesT &&entries) {
            entries_ = std::forward<EntriesT>(entries);
        }

        inline T &at(const SizeType index)
        {
            assert(index < entries_.size());
            return entries_[index];
        }

        inline const T &at(const SizeType index) const
        {
            assert(index < entries_.size());
            return entries_[index];
        }

        void set_rows(SizeType rows) {
            rows_ = rows;
        }

        void set_cols(SizeType cols) {
            cols_ = cols;
        }

        void resize(const SizeType rows, const SizeType cols)
        {
            this->rows_ = rows;
            this->cols_ = cols;
            this->entries_.resize(rows_ * cols_);
        }

        void resize(const Size &s)
        {
            this->rows_ = s.get(0);
            this->cols_ = s.get(1);
            this->entries_.resize(rows_ * cols_);
        }

        T* ptr()
        {
            return &entries_[0];
        }

        const T* ptr() const
        {
            return &entries_[0];
        }

        void identity(const SizeType rows, const SizeType cols, const Scalar &diag = 1.0)
        {
            using std::move;
            using std::min;
            using std::fill;

            resize(rows, cols);
            fill(entries_.begin(), entries_.end(), T(0));

            const SizeType n = min(rows, cols);
            for (SizeType i = 0; i < n; ++i) {
                set(i, i, diag);
            }
        }

        void values(const SizeType rows, const SizeType cols, T &value)
        {
            using std::move;
            using std::min;
            using std::fill;

            resize(rows, cols);
            fill(entries_.begin(), entries_.end(), T(value));
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR MatrixBase ///////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        //locks (No Op)
        inline void read_lock() override {}
        inline void write_lock(WriteMode) override {}
        inline void read_unlock() override {}
        inline void write_unlock(WriteMode) override {}

        //basic mutators
        inline void set(const SizeType &i, const SizeType &j, const T &value) override
        {
            assert(i < rows());
            assert(j < cols());
            entries_[ i + rows() * j ] = value;
        }

        void set(const T &val) override
        {
            std::fill(std::begin(entries_), std::end(entries_), val);
        }

        template<class Index, class Values>
        void set_matrix_aux(
            const Index &rows,
            const Index &cols,
            const Values &vals)
        {
            const SizeType nr = rows.size();
            const SizeType nc = cols.size();
            assert(nr * nc == vals.size());

            SizeType idx = 0;
            for(SizeType i = 0; i < nr; ++i) {
                for(SizeType j = 0; j < nc; ++j) {
                    set(rows[i], rows[j], vals[idx++]);
                }
            }
        }

        //FIXME use span instead
        void set_matrix(
            const std::vector<SizeType> &rows,
            const std::vector<SizeType> &cols,
            const std::vector<T> &vals)
        {
            set_matrix_aux(rows, cols, vals);
        }

        inline void add(const SizeType &i, const SizeType &j, const T &value) override
        {
            assert(i < rows());
            assert(j < cols());
            entries_[ i + rows() * j ] += value;
        }

        //print function
        inline void describe() const override
        {
            auto &os = std::cout;
            describe(os);        }

        inline void describe(std::ostream &os) const
        {

            for (SizeType i=0; i<rows(); ++i)
            {
                for(SizeType j=0; j<cols(); ++j)
                {
                    os << get(i, j) << ' ';
                }
                os << '\n';
            }
        }

        //utility functions
        inline bool empty() const override
        {
            return entries_.empty();
        }

        inline void clear() override
        {
            entries_.clear();
            rows_ = 0;
            cols_ = 0;
        }

        inline Size size() const override
        {
            return {rows_, cols_};
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR DenseMatrix //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline T get(const SizeType &i, const SizeType &j) const override
        {
            assert(i < rows());
            assert(j < cols());
            return at(i + rows() * j);
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR BLAS1Tensor //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////


        ///<T>SWAP - swap x and y
        inline void swap(BlasDenseMatrix &x) override
        {
            std::swap(rows_, x.rows_);
            std::swap(cols_, x.cols_);
            std::swap(entries_, x.entries_);
        }

        ///<T>SCAL - x = a*x
        inline void scale(const T &a) override
        {
            BLASAlgorithms<T>::scal(n_elements(), a, ptr(), 1);
        }

        ///<T>COPY - copy x into y (this)
        inline void copy(const BlasDenseMatrix &x) override
        {
            entries_.resize(x.n_elements());
            rows_ = x.rows_;
            cols_ = x.cols_;
            BLASAlgorithms<T>::copy(x.n_elements(), x.ptr(), 1, ptr(), 1);
        }

        ///<T>AXPY - y = a*x + y
        inline void axpy(const T &a, const BlasDenseMatrix &x) override
        {
            assert(n_elements() == x.n_elements());
            BLASAlgorithms<T>::axpy(n_elements(), a, x.ptr(), 1, ptr(), 1);
        }

        ///<T>DOT - dot product
        inline T dot(const BlasDenseMatrix &other) const override
        {
            assert(n_elements() == other.n_elements());
            return BLASAlgorithms<T>::ddot(n_elements(), ptr(), 1, other.ptr(), 1);
        }

      

        ///<T>AMAX - index of max abs value
        inline SizeType amax() const override
        {
            return BLASAlgorithms<T>::amax(n_elements(), ptr(), 1);
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR BLAS2Matrix //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        void gemv(
            const bool transpose_A,
            const T &alpha,
            const BlasVector &x,
            const T &beta,
            BlasVector &y) const override
        {
            const char t_A_flag = transpose_A ? 'T' : 'N';

            const SizeType m = this->rows();
            const SizeType n = this->cols();

            const SizeType y_size = transpose_A ? this->cols() : this->rows();

            const SizeType lda = m;

            if(y.empty() || y_size != y.size()) {
                y.resize(y_size);
                std::fill(y.begin(), y.end(), 0);
            } else if(approxeq(alpha, 0.0)) {
                std::fill(y.begin(), y.end(), 0);
            }

            BLASAlgorithms<T>::gemv(
                   t_A_flag,               //0
                   m,                      //1
                   n,                      //2
                   alpha,            //3
                   ptr(),               //4
                   lda,                    //5
                   x.ptr(),              //6
                   1,                    //7
                   beta,                   //8
                   y.ptr(),                  //9
                   1
            );                   //10

        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR BLAS3Matrix //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        void gemm(
            const bool transpose_A,
            const T alpha,
            const bool transpose_B,
            const BlasDenseMatrix &B,
            const T beta,
            BlasDenseMatrix &C) const override
        {
            const SizeType k = transpose_A ? this->rows() : this->cols();
            assert(k == (transpose_B ? B.cols() : B.rows()));

            const SizeType m = transpose_A ? this->cols() : this->rows();
            const SizeType n = transpose_B ? B.rows() : B.cols();

            const char t_A_flag = transpose_A ? 'T' : 'N';
            const char t_B_flag = transpose_B ? 'T' : 'N';

            if (C.empty() || approxeq(beta, 0.0)) {
                C.resize(m, n);
                std::fill(C.entries().begin(), C.entries().end(), 0);
            }

            BLASAlgorithms<T>::gemm(
                t_A_flag, 
                t_B_flag, 
                m, 
                n, 
                k, 
                alpha, this->ptr(), 
                transpose_A ? k : m,
                B.ptr(),
                transpose_B ? n : k,
                beta,
                C.ptr(),
                m
            );
        }

        inline void transpose(BlasDenseMatrix &C) const override {
            bool is_squared = rows_ == cols_;

            if(this == &C) {
                if(is_squared) {
                    //in place
                    for(SizeType i = 0; i < rows_; ++i) {
                        for(SizeType j = i + 1; j < rows_; ++j) {
                            C.set(i, j, get(j, i));
                        }
                    }
                } else {
                    BlasDenseMatrix temp;
                    transpose(temp);
                    C = temp;
                }

            } else {
                C.resize(cols_, rows_);

                for(SizeType i = 0; i < rows_; ++i) {
                    for(SizeType j = 0; j < cols_; ++j) {
                        C.set(i, j, get(j, i));
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Constructible //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline void identity(const Size &s, const Scalar &diag = 1.0) override
        {
            assert(s.dims() == 2);
            identity(s.get(0), s.get(1), diag);
        }

        inline void values(const Size &s, const T &val) override 
        {
            resize(s.get(0), s.get(1));
            set(val);
        }

        bool equals(const BlasDenseMatrix &other, const T tol = 0.0) const override
        {
            if(other.rows() != rows()) return false;
            if(other.cols() != cols()) return false;

            const SizeType n = entries_.size();
            assert(n == other.entries_.size());

            for(SizeType i = 0; i < n; ++i) {
                if(std::abs(entries_[i] - other.entries_[i]) > tol) {
                    return false;
                }
            }

            return true;
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR ElementWiseOperand //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////


        inline void e_mul(const BlasDenseMatrix &other) override
        {
            const SizeType n = entries_.size();
            assert(n == other.entries_.size());

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] *= other.entries_[i];
            }  
        }

        inline void e_min(const BlasDenseMatrix &other) override
        {
            const SizeType n = entries_.size();
            assert(n == other.entries_.size());

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] = std::min(other.entries_[i], entries_[i]);
            }  
        }

        inline void e_max(const BlasDenseMatrix &other) override
        {
            const SizeType n = entries_.size();
            assert(n == other.entries_.size());

            for(SizeType i = 0; i < n; ++i) {
                entries_[i] = std::max(other.entries_[i], entries_[i]);
            }  
        }

        ////////////// OVERRIDES FOR Normed //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////

        inline T norm_infty() const override
        {
            using std::max_element;
            if(entries_.empty()) return 0.;
            return *max_element(entries_.begin(), entries_.end());
        }

        ///<T>NRM2 - Euclidean norm
        inline T norm2() const override
        {
            return BLASAlgorithms<T>::nrm2(n_elements(), ptr(), 1);
        }

        ///<T>ASUM - sum of absolute values
        inline T norm1() const override
        {
            return BLASAlgorithms<T>::asum(n_elements(), ptr(), 1);
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Transformable //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////


        void transform(const Sqrt &op) override
        {
            aux_transform(op);
        }

        void transform(const Pow2 &op) override
        {
            aux_transform(op);
        }

        void transform(const Log &op) override
        {
            aux_transform(op);
        }

        void transform(const Exp &op) override
        {
            aux_transform(op);
        }

        void transform(const Cos &op) override
        {
            aux_transform(op);
        }

        void transform(const Sin &op) override
        {
            aux_transform(op);
        }

        void transform(const Abs &op) override
        {
            aux_transform(op);
        }

        void transform(const Pow &op) override
        {
            aux_transform(op);
        }

        void transform(const Reciprocal<T> &op) override
        {
            aux_transform(op);
        }

        void transform(const Minus &) override
        {
            for(auto &e : entries_) {
                e = -e;
            }
        }

        template<class Op>
        void aux_transform(const Op &op)
        {
            for(auto &e : entries_) {
                e = op.apply(e);
            }
        }

        //FIXME make use of interface
        inline void build_diag(BlasVector &v) const 
        {
            auto n = std::min(rows_, cols_);
            v.resize(n);

            for(SizeType i = 0; i < n; ++i) {
                v.set(i, get(i, i));
            }
        }

        inline void diag(const BlasVector &v) 
        {
            auto n = v.size();
            resize(n, n);

            for(SizeType i = 0; i < n; ++i) {
                set(i, i, v.get(i));
            }
        }

        inline void diag(const BlasDenseMatrix &m) 
        {
            auto n = std::min(m.rows_, m.cols_);
            resize(n, n);

            for(SizeType i = 0; i < n; ++i) {
                set(i, i, m.get(i, i));
            }
        }

        inline void shift_diag(const T &val) {
            auto n = std::min(rows_, cols_);

            for(SizeType i = 0; i < n; ++i) {
                ref(i, i) += val;
            }
        }

        inline void select(
            const BlasIndexSet &row_index, 
            const BlasIndexSet &col_index, 
            BlasDenseMatrix &result) const override
        {
            const SizeType r = row_index.size();

            if(!col_index.empty()) {
                const SizeType c = col_index.size();

                result.resize(r, c);

                for(SizeType i = 0; i < r; ++i) {
                    for(SizeType j = 0; j < c; ++j) {
                        result.set(i, j, get(row_index.get(i), col_index.get(j)));
                    }
                }
            } else {
                const SizeType c = cols();

                result.resize(r, c);

                for(SizeType i = 0; i < r; ++i) {
                    for(SizeType j = 0; j < c; ++j) {
                        result.set(i, j, get(row_index.get(i), j));
                    }
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////////
        ////////////// OVERRIDES FOR Reducible //////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////


        inline T reduce(const Min &op) const override
        {
            return aux_reduce(op);
        }

        inline T reduce(const Max &op) const override
        {
            return aux_reduce(op);
        }

        inline T reduce(const Plus &op) const override
        {
            return aux_reduce(op);
        }

        inline SizeType nnz(const Scalar tol = 0.0) const override
        {
            if(entries_.empty()) return 0.0;

            const SizeType n = entries_.size();

            SizeType ret = 0;

            for(SizeType i = 0; i < n; ++i) {
                ret += std::abs(entries_[i]) > tol;
            }

            return ret;
        }

    private:
        Entries entries_;
        SizeType rows_;
        SizeType cols_;

        template<class Op>
        inline T aux_reduce(const Op &op) const
        {
            if(entries_.empty()) return 0.0;

            const SizeType n = entries_.size();

            T ret = entries_[0];

            for(SizeType i = 1; i < n; ++i) {
                ret = op.apply(ret, entries_[i]);
            }

            return ret;
        }

        inline T &ref(const SizeType i, const SizeType j)
        {
            return at(i + rows() * j);
        }
    };

}

#endif //UTOPIA_BLAS_MATRIX_HPP

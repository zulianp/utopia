#ifndef UTOPIA_BLAS_MATRIX_HPP
#define UTOPIA_BLAS_MATRIX_HPP

#include "utopia_Interfaces.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_blas_Traits.hpp"
#include "utopia_BLAS_Operands.hpp"
#include "utopia_blas_Algorithms.hpp"


#include <vector>
#include <memory>
#include <iostream>

namespace utopia {
    template<typename T>
    class BLASDenseMatrix : 
        public DenseMatrix<T, std::size_t>,
        public Tensor<BLASDenseMatrix<T>, 2>,
        public BLAS1Tensor<BLASDenseMatrix<T>>
    {
        typedef std::vector<T> Entries;
        typedef size_t SizeType;

    public:
        BLASDenseMatrix(SizeType rows, SizeType cols) : rows_(rows), cols_(cols)
        {
            entries_.resize(rows_*cols_);
        }

        BLASDenseMatrix(SizeType rows, SizeType cols, T value) : rows_(rows), cols_(cols)
        {
            set_entries(Entries(rows_ * cols_, value));
        }

        ~BLASDenseMatrix() { }

        BLASDenseMatrix() : rows_(0), cols_(0) {}

        BLASDenseMatrix(SizeType rows, SizeType cols, std::initializer_list<T> args) : rows_(rows), cols_(cols)
        {
            using std::copy;

            entries_.resize(rows_*cols_);
            copy(args.begin(), args.end(), entries_.begin());
        }

        BLASDenseMatrix(const Entries& e)
        {
            rows_ = e.size();
            cols_ = 1;
            entries_ = e;

        }
        // inline bool empty() const
        // {
        //     return entries_.empty();
        // }

        SizeType rows() const {
            return rows_;
        }

        SizeType cols() const {
            return cols_;
        }

        SizeType size() const {
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
        // void set(SizeType i, SizeType j, T value) {
        //     assert(i < rows());
        //     assert(j < cols());
        //     entries_[ i + rows() * j ] = value;
        // }

        // T get(SizeType i, SizeType j) const {
        //     assert(i < rows());
        //     assert(j < cols());
        //     return at(i + rows() * j);
        // }

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

        T* ptr()
        {
            return &entries_[0];
        }

        const T* ptr() const
        {
            return &entries_[0];
        }

        void identity(const SizeType rows, const SizeType cols)
        {
            using std::move;
            using std::min;
            using std::fill;

            resize(rows, cols);
            fill(entries_.begin(), entries_.end(), T(0));

            const SizeType n = min(rows, cols);
            for (SizeType i = 0; i < n; ++i)
                set(i, i, T(1.0));
        }

        void values(const SizeType rows, const SizeType cols, T value)
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
        inline void swap(BLASDenseMatrix &x) override
        {
            std::swap(rows_, x.rows_);
            std::swap(cols_, x.cols_);
            std::swap(entries_, x.entries_);
        }

        ///<T>SCAL - x = a*x
        inline void scale(T &a) override
        {
            BLASAlgorithms<T>::scal(size(), a, ptr(), 1);
        }

        ///<T>COPY - copy x into y (this)
        inline void copy(const BLASDenseMatrix &x) override
        {
            entries_.resize(x.size());
            rows_ = x.rows_;
            cols_ = x.cols_;
            BLASAlgorithms<T>::copy(x.size(), x.ptr(), 1, ptr(), 1);
        }

        ///<T>AXPY - y = a*x + y
        inline void axpy(const T &a, const BLASDenseMatrix &x) override
        {
            assert(size() == x.size());
            BLASAlgorithms<T>::axpy(size(), a, x.ptr(), 1, ptr(), 1);
        }

        ///<T>DOT - dot product
        inline T dot(const BLASDenseMatrix &other) const override
        {
            assert(size() == other.size());
            return BLASAlgorithms<T>::ddot(size(), ptr(), 1, other.ptr(), 1);
        }

        ///<T>NRM2 - Euclidean norm
        inline T norm2() const override
        {
            return BLASAlgorithms<T>::nrm2(size(), ptr(), 1);
        }

        ///<T>ASUM - sum of absolute values
        inline T asum() const override
        {
            return BLASAlgorithms<T>::asum(size(), ptr(), 1);
        }

        ///<T>AMAX - index of max abs value
        inline SizeType amax() const override
        {
            return BLASAlgorithms<T>::amax(size(), ptr(), 1);
        }

    private:
        Entries entries_;
        SizeType rows_;
        SizeType cols_;

    };

    inline Wrapper<BLASDenseMatrix<double>, 2> mmake(int rows, int cols, std::initializer_list<double> args) {
        return BLASDenseMatrix<double>(rows, cols, args);
    }

    template<typename T>
    void disp(const Wrapper< BLASDenseMatrix<T>, 2> &w, std::ostream &os)
    {
        w.implemenetation().describe(os);
    }

}

#endif //UTOPIA_BLAS_MATRIX_HPP

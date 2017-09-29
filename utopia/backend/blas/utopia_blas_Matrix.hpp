#ifndef UTOPIA_BLAS_MATRIX_HPP
#define UTOPIA_BLAS_MATRIX_HPP

#include "utopia_Wrapper.hpp"
#include "utopia_blas_Traits.hpp"

#include <vector>
#include <memory>
#include <iostream>

namespace utopia {
    template<typename T>
    class Matrix
    {
        typedef std::vector<T> Entries;
        typedef typename Entries::size_type SizeType;

    public:
        Matrix(SizeType rows, SizeType cols) : rows_(rows), cols_(cols)
        {
            entries_.resize(rows_*cols_);
        }

        Matrix(SizeType rows, SizeType cols, T value) : rows_(rows), cols_(cols)
        {
            set_entries(Entries(rows_ * cols_, value));
        }
		
		~Matrix() { }

        Matrix() : rows_(0), cols_(0) {}

        Matrix(SizeType rows, SizeType cols, std::initializer_list<T> args) : rows_(rows), cols_(cols)
        {
            using std::copy;

            entries_.resize(rows_*cols_);
            copy(args.begin(), args.end(), entries_.begin());
        }

        Matrix(const Entries& e)
        {
            rows_ = e.size();
            cols_ = 1;
            entries_ = e;

        }
        inline bool empty() const
        {
            return entries_.empty();
        }

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
        void set(SizeType i, SizeType j, T value) {
            assert(i < rows());
            assert(j < cols());
            entries_[ i + rows() * j ] = value;
        }

        T get(SizeType i, SizeType j) const {
            assert(i < rows());
            assert(j < cols());
            return at(i + rows() * j);
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

    private:
        Entries entries_;
        SizeType rows_;
        SizeType cols_;

    };

    inline Wrapper<Matrix<double>, 2> mmake(int rows, int cols, std::initializer_list<double> args) {
        return Matrix<double>(rows, cols, args);
    }

    template<typename T>
    void disp(const Wrapper< Matrix<T>, 2> &w, std::ostream &os)
    {
        for (SizeType i=0; i<w.implementation().rows(); ++i)
        {
            for(SizeType j=0; j<w.implementation().cols(); ++j)
            {
                os << w.implementation().get(i, j) << ' ';
            }
            os << '\n';
        }
        //disp(w.implementation().entries().begin(), w.implementation().entries().end(), os);
    }

}

#endif //UTOPIA_BLAS_MATRIX_HPP

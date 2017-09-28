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
        Matrix(SizeType rows, SizeType cols) : rows(rows), cols(cols)
        {
            entries.resize(rows*cols);
        }

        Matrix(SizeType rows, SizeType cols, T value) : rows(rows), cols(cols)
        {
            setEntries(Entries(rows*cols, value));
        }
		
		~Matrix()
		{
		}

        Matrix() : rows(0), cols(0) {}

        Matrix(SizeType rows, SizeType cols, std::initializer_list<T> args) : rows(rows), cols(cols)
        {
            entries.resize(rows*cols);
            using std::copy;
            copy(args.begin(), args.end(), entries.begin());
        }

        Matrix(const Entries& e)
        {
            rows = e.size();
            cols = 1;
            entries = e;

        }

        SizeType getRows() const {
            return rows;
        }

        SizeType getCols() const {
            return cols;
        }

        SizeType size() const {
            return getRows() * getCols();
        }

        const Entries &getEntries() const {
            return entries;
        }

        Entries &getEntries() {
            return entries;
        }

        template<class EntriesT>
        void setEntries(EntriesT &&entries) {
            Matrix::entries = std::forward<EntriesT>(entries);
        }
		
		inline T &at(const SizeType index)
		{
			assert(index < entries.size());
			return entries[index];
		}

		inline const T &at(const SizeType index) const
		{
			assert(index < entries.size());
			return entries[index];
		}
        void set(SizeType i, SizeType j, T value) {
            assert(i < getRows());
            assert(j < getCols());
            entries[ i + getRows() * j ] = value;
        }

        T get(SizeType i, SizeType j) const {
            assert(i < getRows());
            assert(j < getCols());
            return at(i + getRows() * j);
        }

        void setRows(SizeType rows) {
            Matrix::rows = rows;
        }

        void setCols(SizeType cols) {
            Matrix::cols = cols;
        }

        void resize(const SizeType rows, const SizeType cols)
        {
            this->rows = rows;
            this->cols = cols;
            this->entries.resize(rows*cols);
        }

        T* ptr()
        {
            return &entries[0];
        }

        const T* ptr() const
        {
            return &entries[0];
        }

        void identity(const SizeType rows, const SizeType cols)
        {
            using std::move;
            using std::min;
            using std::fill;

            resize(rows, cols);
            fill(entries.begin(), entries.end(), T(0));

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
            fill(entries.begin(), entries.end(), T(value));
        }

    private:
        Entries entries;
        SizeType rows;
        SizeType cols;

    };


    inline Wrapper<Matrix<double>, 2> mmake(int rows, int cols, std::initializer_list<double> args) {
        return Matrix<double>(rows, cols, args);
    }



    template<typename T>
    void disp(const Wrapper< Matrix<T>, 2> &w, std::ostream &os)
    {
        for (SizeType i=0; i<w.implementation().getRows(); ++i)
        {
            for(SizeType j=0; j<w.implementation().getCols(); ++j)
            {
                os << w.implementation().get(i, j) << ' ';
            }
            os << '\n';
        }
        //disp(w.implementation().getEntries().begin(), w.implementation().getEntries().end(), os);
    }

}

#endif //UTOPIA_BLAS_MATRIX_HPP

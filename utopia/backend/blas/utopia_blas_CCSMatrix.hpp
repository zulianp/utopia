
#ifndef UTOPIA_UTOPIA_CCSMATRIX_HPP
#define UTOPIA_UTOPIA_CCSMATRIX_HPP

#include <vector>
#include <map>

namespace utopia {
    template<typename Scalar>
    class CCSMatrix {
    public:
        typedef std::vector<Scalar> Entries;
        typedef int SizeType;
        typedef std::map< std::pair<SizeType, SizeType>, Scalar > MapMatrix;
        static const SizeType INVALID_INDEX = -1;

        template<typename SizeTypeT1, typename SizeTypeT2, typename ScalarT>
        CCSMatrix(const SizeType rows, const SizeType cols,
                  std::initializer_list<SizeTypeT1> colptr,
                  std::initializer_list<SizeTypeT2> rowindex,
                  std::initializer_list<ScalarT> entries)
                : rows_(rows), cols_(cols),
                  colptr_(colptr.size()),
                  rowindex_(rowindex.size()),
                  entries_(entries.size()),
                  _editing(false)
        {
            using std::copy;
            copy(colptr.begin(), colptr.end(), colptr_.begin());
            copy(rowindex.begin(), rowindex.end(), rowindex_.begin());
            copy(entries.begin(), entries.end(), entries_.begin());
        }

        void initialize(const SizeType rows, const SizeType cols, const SizeType nnzXCol)
        {
            using std::fill;
            rows_ = rows;
            cols_ = cols;

            assert(nnzXCol>0);

            colptr_.resize(cols_+1);
            rowindex_.resize(rows_*nnzXCol);
            entries_.resize(rowindex_.size());

            fill(entries_.begin(), entries_.end(), 0);

            colptr_[0] = 0;

            const SizeType col_ptr_size = colptr_.size();
            for(SizeType i = 1; i < col_ptr_size; ++i) {
                colptr_[i] = nnzXCol + colptr_[i-1];
            }

            for(SizeType i = 0; i < rows_; ++i) {
                for(SizeType k = 0; k < nnzXCol; ++k) {
                    rowindex_[i*nnzXCol+k] = k;
                }
            }

            _editing = false;
        }

        template<typename SizeTypeT1, typename SizeTypeT2, typename ScalarT>
        void initialize(const SizeType rows, const SizeType cols,
                        std::initializer_list<SizeTypeT1> colptr,
                        std::initializer_list<SizeTypeT2> rowindex,
                        std::initializer_list<ScalarT> entries)
        {
            using std::copy;
            rows_ = rows;
            cols_ = cols;

            colptr_.resize(rows+1);
            rowindex_.resize(rowindex.size());
            entries_.resize(entries.size());

            copy(colptr.begin(), colptr.end(), colptr_.begin());
            copy(rowindex.begin(), rowindex.end(), rowindex_.begin());
            copy(entries.begin(), entries.end(), entries_.begin());
        }

        CCSMatrix(const SizeType rows, const SizeType cols, const SizeType nEntries)
                : rows_(rows), cols_(cols), colptr_(rows+1, 0), rowindex_(nEntries, 0), entries_(nEntries, 0),  _editing(false)
        {}

        CCSMatrix()
                : rows_(0), cols_(0), _editing(false)
        {}

        void clear() {
            rows_ = 0;
            cols_ = 0;
            colptr_.clear();
            rowindex_.clear();
            entries_.clear();
            _buffer.clear();
            _editing = false;
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

        template<class EntriesT>
        void set_entries(EntriesT &&entries) {
            entries_ = std::forward<EntriesT>(entries);
        }

        inline Scalar &at(const SizeType index)
        {
            assert(!_editing);
            assert(index < entries_.size());
            return entries_[index];
        }

        inline const Scalar &at(const SizeType index) const
        {
            assert(!_editing);
            assert(index < entries_.size());
            return entries_[index];
        }

        void assembly_begin()
        {
            _editing = true;
        }

        /// Triggers a clear();
        void resize(const SizeType rows, const SizeType cols)
        {
            clear();
            this->rows_ = rows;
            this->cols_ = cols;
        }

        void from(const MapMatrix &buffer)
        {
            colptr_.resize(rows()+1);
            fill(colptr_.begin(), colptr_.end(), 0);
            rowindex_.resize(buffer.size());
            entries_.resize(buffer.size());


            for(auto &e : buffer) {
                ++colptr_[e.first.first+1];
            }

            for(SizeType i = 0; i < cols(); ++i) {
                colptr_[i+1] += colptr_[i];
            }

            std::vector<SizeType> col_offsets(cols());
            fill(col_offsets.begin(), col_offsets.end(), 0);

            for(auto &e : buffer) {
                assert(e.first.first >= 0);
                assert(e.first.second >= 0);
                assert(e.first.first  < cols());
                assert(e.first.second < rows());

                const SizeType index = colptr_[e.first.first] + col_offsets[e.first.first];
                rowindex_[index] = e.first.second;
                entries_[index] = e.second;
                ++col_offsets[e.first.first];
            }
        }

        void assembly_end()
        {
            using std::fill;

            if(!_editing) {
                assert(_buffer.empty());
                return;
            }

            if(_buffer.empty())
                return;

            put(_buffer);
            from(_buffer);

            _buffer.clear();
            _editing = false;
        }

        void put(MapMatrix &mat, const bool map_transpose = false)
        {
            for(SizeType c = 0; c != rows(); ++c) {
                for(SizeType k = colptr_[c]; k != colptr_[c+1]; ++k) {
                    const SizeType r = rowindex_[k];
                    if(map_transpose) {
                        mat[std::make_pair(r, c)] = entries_[k];
                    } else {
                        mat[std::make_pair(c, r)] = entries_[k];
                    }
                }
            }
        }


        void set(const SizeType i, const SizeType j, const Scalar value) {
            assert(i < rows());
            assert(j < cols());
            assert(_editing);

            if(value == 0.0) return;

            const SizeType index = find_efficient(i, j);
            assert(find(i, j) == index);
            if(index != INVALID_INDEX) {
                entries_[index] = value;
                return;
            }

            _buffer[std::make_pair(j, i)] = value;
        }

        inline SizeType find_efficient(const SizeType i, const SizeType j) const
        {
            using std::lower_bound;
            using std::distance;

            if(colptr_.empty()) return INVALID_INDEX;

            if(colptr_[j+1] >= rowindex_.size()) {
                return INVALID_INDEX;
            }

            auto low = lower_bound(rowindex_.begin() + colptr_[j],  rowindex_.begin() + colptr_[j+1], j);
            return (*low == i)? SizeType(distance(rowindex_.begin(), low)) : INVALID_INDEX;
        }


        SizeType find(const SizeType i, const SizeType j)
        {
            if(colptr_.empty()) return INVALID_INDEX;

            for(SizeType k = colptr_[j]; k != colptr_[j+1]; ++k) {
                if(i == rowindex_[k]) {
                    return k;
                }
            }

            return INVALID_INDEX;
        }

        Scalar get(const SizeType i, const SizeType j) const {
            assert(!_editing);
            for(SizeType k = colptr_[j]; k != colptr_[j+1]; ++k) {
                if(i == rowindex_[k]) {
                    return entries_[k];
                }
            }

            return 0;
        }

        void setRows(SizeType rows) {
            rows_ = rows;
        }

        void setCols(SizeType cols) {
            cols_ = cols;
        }

        inline const std::vector<SizeType> &colptr() const {
            return colptr_;
        }

        inline const std::vector<SizeType> &rowindex() const {
            return rowindex_;
        }



    private:
        SizeType rows_;
        SizeType cols_;
        std::vector<SizeType> colptr_;
        std::vector<SizeType> rowindex_;
        std::vector<Scalar> entries_;

        MapMatrix _buffer;
        bool _editing;
    };


    template<typename T>
    void disp(const Wrapper< CCSMatrix<T>, 2> &w, std::ostream &os)
    {
        os << "-------------------------------------------------------------------------\n";
        os << w.implementation().rows() << ", " << w.implementation().cols() << "\n";
        os << "colptr:\n";
        disp(w.implementation().colptr().begin(), w.implementation().colptr().end(), os);
        os << "\nrowindex:\n";
        disp(w.implementation().rowindex().begin(), w.implementation().rowindex().end(), os);
        os << "\nentries:\n";
        disp(w.implementation().entries().begin(), w.implementation().entries().end(), os);
        os << "-------------------------------------------------------------------------\n";

    }
}


#endif //UTOPIA_UTOPIA_CCSMATRIX_HPP

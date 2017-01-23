
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
                : _rows(rows), _cols(cols),
                  _colptr(colptr.size()),
                  _rowindex(rowindex.size()),
                  _entries(entries.size()),
                  _editing(false)
        {
            using std::copy;
            copy(colptr.begin(), colptr.end(), _colptr.begin());
            copy(rowindex.begin(), rowindex.end(), _rowindex.begin());
            copy(entries.begin(), entries.end(), _entries.begin());
        }

        void initialize(const SizeType rows, const SizeType cols, const SizeType nnzXCol)
        {
            using std::fill;
            _rows = rows;
            _cols = cols;

            assert(nnzXCol>0);

            _colptr.resize(_cols+1);
            _rowindex.resize(_rows*nnzXCol);
            _entries.resize(_rowindex.size());

            fill(_entries.begin(), _entries.end(), 0);

            _colptr[0] = 0;
            for(SizeType i = 1; i < _colptr.size(); ++i) {
                _colptr[i] = nnzXCol + _colptr[i-1];
            }

            for(SizeType i = 0; i < _rows; ++i) {
                for(SizeType k = 0; k < nnzXCol; ++k) {
                    _rowindex[i*nnzXCol+k] = k;
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
            _rows = rows;
            _cols = cols;

            _colptr.resize(rows+1);
            _rowindex.resize(rowindex.size());
            _entries.resize(entries.size());

            copy(colptr.begin(), colptr.end(), _colptr.begin());
            copy(rowindex.begin(), rowindex.end(), _rowindex.begin());
            copy(entries.begin(), entries.end(), _entries.begin());
        }

        CCSMatrix(const SizeType rows, const SizeType cols, const SizeType nEntries)
                : _rows(rows), _cols(cols), _colptr(rows+1, 0), _rowindex(nEntries, 0), _entries(nEntries, 0),  _editing(false)
        {}

        CCSMatrix()
                : _rows(0), _cols(0), _editing(false)
        {}

        void clear() {
            _rows = 0;
            _cols = 0;
            _colptr.clear();
            _rowindex.clear();
            _entries.clear();
            _buffer.clear();
            _editing = false;
        }

        SizeType getRows() const {
            return _rows;
        }

        SizeType getCols() const {
            return _cols;
        }

        SizeType size() const {
            return getRows() * getCols();
        }

        const Entries &getEntries() const {
            return _entries;
        }

        template<class EntriesT>
        void setEntries(EntriesT &&entries) {
            _entries = std::forward<EntriesT>(entries);
        }

        inline Scalar &at(const SizeType index)
        {
            assert(!_editing);
            assert(index < _entries.size());
            return _entries[index];
        }

        inline const Scalar &at(const SizeType index) const
        {
            assert(!_editing);
            assert(index < _entries.size());
            return _entries[index];
        }

        void assemblyBegin()
        {
            _editing = true;
        }

        /// Triggers a clear();
        void resize(const SizeType rows, const SizeType cols)
        {
            clear();
            this->_rows = rows;
            this->_cols = cols;
        }

        void from(const MapMatrix &buffer)
        {
            _colptr.resize(getRows()+1);
            fill(_colptr.begin(), _colptr.end(), 0);
            _rowindex.resize(buffer.size());
            _entries.resize(buffer.size());


            for(auto &e : buffer) {
                ++_colptr[e.first.first+1];
            }

            for(SizeType i = 0; i < getCols(); ++i) {
                _colptr[i+1] += _colptr[i];
            }

            std::vector<SizeType> col_offsets(getCols());
            fill(col_offsets.begin(), col_offsets.end(), 0);

            for(auto &e : buffer) {
                assert(e.first.first >= 0);
                assert(e.first.second >= 0);
                assert(e.first.first  < getCols());
                assert(e.first.second < getRows());

                const SizeType index = _colptr[e.first.first] + col_offsets[e.first.first];
                _rowindex[index] = e.first.second;
                _entries[index] = e.second;
                ++col_offsets[e.first.first];
            }
        }

        void assemblyEnd()
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
            for(SizeType c = 0; c != getRows(); ++c) {
                for(SizeType k = _colptr[c]; k != _colptr[c+1]; ++k) {
                    const SizeType r = _rowindex[k];
                    if(map_transpose) {
                        mat[std::make_pair(r, c)] = _entries[k];
                    } else {
                        mat[std::make_pair(c, r)] = _entries[k];
                    }
                }
            }
        }


        void set(const SizeType i, const SizeType j, const Scalar value) {
            assert(i < getRows());
            assert(j < getCols());
            assert(_editing);

            if(value == 0.0) return;

            const SizeType index = find_efficient(i, j);
            assert(find(i, j) == index);
            if(index != INVALID_INDEX) {
                _entries[index] = value;
                return;
            }

            _buffer[std::make_pair(j, i)] = value;
        }

        inline SizeType find_efficient(const SizeType i, const SizeType j) const
        {
            using std::lower_bound;
            using std::distance;

            if(_colptr.empty()) return INVALID_INDEX;
            auto low = lower_bound(_rowindex.begin() + _colptr[j],  _rowindex.begin() + _colptr[j+1], j);
            return (*low == i)? SizeType(distance(_rowindex.begin(), low)) : INVALID_INDEX;
        }


        SizeType find(const SizeType i, const SizeType j)
        {
            if(_colptr.empty()) return INVALID_INDEX;

            for(SizeType k = _colptr[j]; k != _colptr[j+1]; ++k) {
                if(i == _rowindex[k]) {
                    return k;
                }
            }

            return INVALID_INDEX;
        }

        Scalar get(const SizeType i, const SizeType j) const {
            assert(!_editing);
            for(SizeType k = _colptr[j]; k != _colptr[j+1]; ++k) {
                if(i == _rowindex[k]) {
                    return _entries[k];
                }
            }

            return 0;
        }

        void setRows(SizeType rows) {
            _rows = rows;
        }

        void setCols(SizeType cols) {
            _cols = cols;
        }

        inline const std::vector<SizeType> &colptr() const {
            return _colptr;
        }

        inline const std::vector<SizeType> &rowindex() const {
            return _rowindex;
        }



    private:
        SizeType _rows;
        SizeType _cols;
        std::vector<SizeType> _colptr;
        std::vector<SizeType> _rowindex;
        std::vector<Scalar> _entries;

        MapMatrix _buffer;
        bool _editing;
    };


    template<typename T>
    void disp(const Wrapper< CCSMatrix<T>, 2> &w, std::ostream &os)
    {
        os << "-------------------------------------------------------------------------\n";
        os << w.implementation().getRows() << ", " << w.implementation().getCols() << "\n";
        os << "colptr:\n";
        disp(w.implementation().colptr().begin(), w.implementation().colptr().end(), os);
        os << "\nrowindex:\n";
        disp(w.implementation().rowindex().begin(), w.implementation().rowindex().end(), os);
        os << "\nentries:\n";
        disp(w.implementation().getEntries().begin(), w.implementation().getEntries().end(), os);
        os << "-------------------------------------------------------------------------\n";

    }
}


#endif //UTOPIA_UTOPIA_CCSMATRIX_HPP

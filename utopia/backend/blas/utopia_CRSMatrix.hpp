
#ifndef UTOPIA_UTOPIA_CRSMATRIX_HPP
#define UTOPIA_UTOPIA_CRSMATRIX_HPP

#include <vector>
#include <map>

namespace utopia {
    template<typename Scalar>
    class CRSMatrix {
    public:
        typedef std::vector<Scalar> Entries;
        typedef int SizeType;
        typedef std::map< std::pair<SizeType, SizeType>, Scalar > MapMatrix;
        static const SizeType INVALID_INDEX = -1;

        template<typename SizeTypeT1, typename SizeTypeT2, typename ScalarT>
        CRSMatrix(const SizeType rows, const SizeType cols,
                  std::initializer_list<SizeTypeT1> rowptr,
                  std::initializer_list<SizeTypeT2> colindex,
                  std::initializer_list<ScalarT> entries)
                : _rows(rows), _cols(cols),
                  _rowptr(rowptr.size()),
                  _colindex(colindex.size()),
                  _entries(entries.size()),
                  _editing(false)
        {
            using std::copy;
            copy(rowptr.begin(), rowptr.end(), _rowptr.begin());
            copy(colindex.begin(), colindex.end(), _colindex.begin());
            copy(entries.begin(), entries.end(), _entries.begin());
        }

        template<typename SizeTypeT1, typename SizeTypeT2, typename ScalarT>
        void initialize(const SizeType rows, const SizeType cols,
                        std::initializer_list<SizeTypeT1> rowptr,
                        std::initializer_list<SizeTypeT2> colindex,
                        std::initializer_list<ScalarT> entries)
        {
            using std::copy;
            _rows = rows;
            _cols = cols;

            _rowptr.resize(rows+1);
            _colindex.resize(colindex.size());
            _entries.resize(entries.size());

            copy(rowptr.begin(), rowptr.end(), _rowptr.begin());
            copy(colindex.begin(), colindex.end(), _colindex.begin());
            copy(entries.begin(), entries.end(), _entries.begin());
        }

        CRSMatrix(const SizeType rows, const SizeType cols, const SizeType nEntries)
        : _rows(rows), _cols(cols), _rowptr(rows+1, 0), _colindex(nEntries, 0), _entries(nEntries, 0),  _editing(false)
        {}

        CRSMatrix()
        : _rows(0), _cols(0), _editing(false)
        {}


        void initialize(const SizeType rows, const SizeType cols, const SizeType nnzXRow)
        {
            using std::fill;
            _rows = rows;
            _cols = cols;

            assert(nnzXRow>0);

            _rowptr.resize(_rows+1);
            _colindex.resize(_cols*nnzXRow);
            _entries.resize(_colindex.size());

            fill(_entries.begin(), _entries.end(), 0);

            _rowptr[0] = 0;
            for(SizeType i = 1; i < _rowptr.size(); ++i) {
                _rowptr[i] = nnzXRow + _rowptr[i-1];
            }

            for(SizeType i = 0; i < _cols; ++i) {
                for(SizeType k = 0; k < nnzXRow; ++k) {
                    _colindex[i*nnzXRow+k] = k;
                }
            }

            _editing = false;
        }

        void clear() {
            _rows = 0;
            _cols = 0;
            _rowptr.clear();
            _colindex.clear();
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

        ///! Triggers a clear
        void resize(const SizeType rows, const SizeType cols)
        {
            clear();
            _rows = rows;
            _cols = cols;
        }

        const Entries &getEntries() const {
            return _entries;
        }

        Entries &getEntries() {
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

            _rowptr.resize(getRows()+1);
           fill(_rowptr.begin(), _rowptr.end(), 0);
            _colindex.resize(_buffer.size());
            _entries.resize(_buffer.size());


            for(auto &e : _buffer) {
                ++_rowptr[e.first.first+1];
            }

            for(SizeType i = 0; i < getRows(); ++i) {
                _rowptr[i+1] += _rowptr[i];
            }

            std::vector<SizeType> row_offsets(getRows());
            fill(row_offsets.begin(), row_offsets.end(), 0);

            for(auto &e : _buffer) {
                assert(e.first.first >= 0);
                assert(e.first.second >= 0);
                assert(e.first.first < getRows());
                assert(e.first.second < getCols());

                const SizeType index = _rowptr[e.first.first] + row_offsets[e.first.first];
                _colindex[index] = e.first.second;
                _entries[index] = e.second;
                ++row_offsets[e.first.first];
            }

            _buffer.clear();
            _editing = false;
        }

        void put(MapMatrix &mat, const bool map_transpose = false) const
        {
            if(_rowptr.empty()) return;

            for(SizeType r = 0; r < getRows(); ++r) {
                for(SizeType k = _rowptr[r]; k != _rowptr[r+1]; ++k) {
                    const SizeType c = _colindex[k];
                    if(map_transpose) {
                        mat[std::make_pair(c, r)] = _entries[k];
                    } else {
                        mat[std::make_pair(r, c)] = _entries[k];
                    }
                }
            }
        }


        void set(const SizeType i, const SizeType j, const Scalar value) {
            assert(i < getRows());
            assert(j < getCols());
            assert(_editing);

            if(value == 0.0) return;

//            const SizeType index = find(i, j);
            const SizeType index = find_efficient(i, j);
            assert(find(i, j) == index);

            if(index != INVALID_INDEX) {
                assert(index < _entries.size());
                _entries[index] = value;
                return;
            }

            _buffer[std::make_pair(i, j)] = value;
        }

        inline SizeType find_efficient(const SizeType i, const SizeType j) const
        {
            using std::lower_bound;
            using std::distance;

            if(_rowptr.empty()) return INVALID_INDEX;

            assert(i + 1 < _rowptr.size());
            
            if(_rowptr[i + 1] >= _colindex.size()) {
                return INVALID_INDEX;
            }

            auto low = lower_bound(_colindex.begin() + _rowptr[i],  _colindex.begin() + _rowptr[i + 1], j);
            return (*low == j)? SizeType(distance(_colindex.begin(), low)) : INVALID_INDEX;
        }

        inline SizeType find(const SizeType i, const SizeType j) const {
            if(_rowptr.empty()) return INVALID_INDEX;
            for(SizeType k = _rowptr[i]; k != _rowptr[i+1]; ++k) {
                if(j == _colindex[k]) {
                    return k;
                }
            }

            return INVALID_INDEX;
        }

        Scalar get(const SizeType i, const SizeType j) const {
            assert(!_editing);
            for(SizeType k = _rowptr[i]; k != _rowptr[i+1]; ++k) {
                if(j == _colindex[k]) {
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

        inline const std::vector<SizeType> &rowptr() const {
            return _rowptr;
        }

        inline const std::vector<SizeType> &colindex() const {
            return _colindex;
        }



    private:
        SizeType _rows;
        SizeType _cols;
        std::vector<SizeType> _rowptr;
        std::vector<SizeType> _colindex;
        std::vector<Scalar> _entries;

        MapMatrix _buffer;
        bool _editing;
    };


    template<typename T>
    void disp(const Wrapper< CRSMatrix<T>, 2> &w, std::ostream &os)
    {
        os << "-------------------------------------------------------------------------\n";
        os << w.implementation().getRows() << ", " << w.implementation().getCols() << "\n";
        os << "rowptr:\n";
        disp(w.implementation().rowptr().begin(), w.implementation().rowptr().end(), os);
        os << "\ncolindex:\n";
        disp(w.implementation().colindex().begin(), w.implementation().colindex().end(), os);
        os << "\nentries:\n";
        disp(w.implementation().getEntries().begin(), w.implementation().getEntries().end(), os);
        os << "-------------------------------------------------------------------------\n";

    }
}

#endif //UTOPIA_UTOPIA_CRSMATRIX_HPP


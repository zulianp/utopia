
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
                : rows_(rows), cols_(cols),
                  rowptr_(rowptr.size()),
                  colindex_(colindex.size()),
                  entries_(entries.size()),
                  editing_(false)
        {
            using std::copy;
            copy(rowptr.begin(), rowptr.end(), rowptr_.begin());
            copy(colindex.begin(), colindex.end(), colindex_.begin());
            copy(entries.begin(), entries.end(), entries_.begin());
        }

        template<typename SizeTypeT1, typename SizeTypeT2, typename ScalarT>
        void initialize(const SizeType rows, const SizeType cols,
                        std::initializer_list<SizeTypeT1> rowptr,
                        std::initializer_list<SizeTypeT2> colindex,
                        std::initializer_list<ScalarT> entries)
        {
            using std::copy;
            rows_ = rows;
            cols_ = cols;

            rowptr_.resize(rows+1);
            colindex_.resize(colindex.size());
            entries_.resize(entries.size());

            copy(rowptr.begin(), rowptr.end(), rowptr_.begin());
            copy(colindex.begin(), colindex.end(), colindex_.begin());
            copy(entries.begin(), entries.end(), entries_.begin());
        }

        CRSMatrix(const SizeType rows, const SizeType cols, const SizeType nEntries)
        : rows_(rows), cols_(cols), rowptr_(rows+1, 0), colindex_(nEntries, 0), entries_(nEntries, 0),  editing_(false)
        {}

        CRSMatrix()
        : rows_(0), cols_(0), editing_(false)
        {}


        void initialize(const SizeType rows, const SizeType cols, const SizeType nnzXRow)
        {
            using std::fill;
            rows_ = rows;
            cols_ = cols;

            assert(nnzXRow>0);

            rowptr_.resize(rows_+1);
            colindex_.resize(cols_*nnzXRow);
            entries_.resize(colindex_.size());

            fill(entries_.begin(), entries_.end(), 0);

            rowptr_[0] = 0;
            for(SizeType i = 1; i < rowptr_.size(); ++i) {
                rowptr_[i] = nnzXRow + rowptr_[i-1];
            }

            for(SizeType i = 0; i < cols_; ++i) {
                for(SizeType k = 0; k < nnzXRow; ++k) {
                    colindex_[i*nnzXRow+k] = k;
                }
            }

            editing_ = false;
        }

        void clear() {
            rows_ = 0;
            cols_ = 0;
            rowptr_.clear();
            colindex_.clear();
            entries_.clear();
            buffer_.clear();
            editing_ = false;
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

        ///! Triggers a clear
        void resize(const SizeType rows, const SizeType cols)
        {
            clear();
            rows_ = rows;
            cols_ = cols;
        }

        const Entries &entries() const {
            return entries_;
        }

        Entries &entries() {
            return entries_;
        }

        template<class EntriesT>
        void setEntries(EntriesT &&entries) {
            entries_ = std::forward<EntriesT>(entries);
        }

        inline Scalar &at(const SizeType index)
        {
            assert(!editing_);
            assert(index < entries_.size());
            return entries_[index];
        }

        inline const Scalar &at(const SizeType index) const
        {
            assert(!editing_);
            assert(index < entries_.size());
            return entries_[index];
        }

        void assembly_begin()
        {
            editing_ = true;
        }

        void assembly_end()
        {
            using std::fill;

            if(!editing_) {
                assert(buffer_.empty());
                return;
            }

            if(buffer_.empty())
                return;

            put(buffer_);

            rowptr_.resize(rows()+1);
           fill(rowptr_.begin(), rowptr_.end(), 0);
            colindex_.resize(buffer_.size());
            entries_.resize(buffer_.size());


            for(auto &e : buffer_) {
                ++rowptr_[e.first.first+1];
            }

            for(SizeType i = 0; i < rows(); ++i) {
                rowptr_[i+1] += rowptr_[i];
            }

            std::vector<SizeType> row_offsets(rows());
            fill(row_offsets.begin(), row_offsets.end(), 0);

            for(auto &e : buffer_) {
                assert(e.first.first >= 0);
                assert(e.first.second >= 0);
                assert(e.first.first < rows());
                assert(e.first.second < cols());

                const SizeType index = rowptr_[e.first.first] + row_offsets[e.first.first];
                colindex_[index] = e.first.second;
                entries_[index] = e.second;
                ++row_offsets[e.first.first];
            }

            buffer_.clear();
            editing_ = false;
        }

        void put(MapMatrix &mat, const bool map_transpose = false) const
        {
            if(rowptr_.empty()) return;

            for(SizeType r = 0; r < rows(); ++r) {
                for(SizeType k = rowptr_[r]; k != rowptr_[r+1]; ++k) {
                    const SizeType c = colindex_[k];
                    if(map_transpose) {
                        mat[std::make_pair(c, r)] = entries_[k];
                    } else {
                        mat[std::make_pair(r, c)] = entries_[k];
                    }
                }
            }
        }


        void set(const SizeType i, const SizeType j, const Scalar value) {
            assert(i < rows());
            assert(j < cols());
            assert(editing_);

            if(value == 0.0) return;

//            const SizeType index = find(i, j);
            const SizeType index = find_efficient(i, j);
            assert(find(i, j) == index);

            if(index != INVALID_INDEX) {
                assert(index < entries_.size());
                entries_[index] = value;
                return;
            }

            buffer_[std::make_pair(i, j)] = value;
        }

        inline SizeType find_efficient(const SizeType i, const SizeType j) const
        {
            using std::lower_bound;
            using std::distance;

            if(rowptr_.empty()) return INVALID_INDEX;

            assert(i + 1 < rowptr_.size());
            
            if(rowptr_[i + 1] >= colindex_.size()) {
                return INVALID_INDEX;
            }

            auto low = lower_bound(colindex_.begin() + rowptr_[i],  colindex_.begin() + rowptr_[i + 1], j);
            return (*low == j)? SizeType(distance(colindex_.begin(), low)) : INVALID_INDEX;
        }

        inline SizeType find(const SizeType i, const SizeType j) const {
            if(rowptr_.empty()) return INVALID_INDEX;
            for(SizeType k = rowptr_[i]; k != rowptr_[i+1]; ++k) {
                if(j == colindex_[k]) {
                    return k;
                }
            }

            return INVALID_INDEX;
        }

        Scalar get(const SizeType i, const SizeType j) const {
            assert(!editing_);
            for(SizeType k = rowptr_[i]; k != rowptr_[i+1]; ++k) {
                if(j == colindex_[k]) {
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

        inline const std::vector<SizeType> &rowptr() const {
            return rowptr_;
        }

        inline const std::vector<SizeType> &colindex() const {
            return colindex_;
        }



    private:
        SizeType rows_;
        SizeType cols_;
        std::vector<SizeType> rowptr_;
        std::vector<SizeType> colindex_;
        std::vector<Scalar> entries_;

        MapMatrix buffer_;
        bool editing_;
    };


    template<typename T>
    void disp(const Wrapper< CRSMatrix<T>, 2> &w, std::ostream &os)
    {
        os << "-------------------------------------------------------------------------\n";
        os << w.implementation().rows() << ", " << w.implementation().cols() << "\n";
        os << "rowptr:\n";
        disp(w.implementation().rowptr().begin(), w.implementation().rowptr().end(), os);
        os << "\ncolindex:\n";
        disp(w.implementation().colindex().begin(), w.implementation().colindex().end(), os);
        os << "\nentries:\n";
        disp(w.implementation().entries().begin(), w.implementation().entries().end(), os);
        os << "-------------------------------------------------------------------------\n";

    }
}

#endif //UTOPIA_UTOPIA_CRSMATRIX_HPP


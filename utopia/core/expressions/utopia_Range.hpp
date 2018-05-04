#ifndef UTOPIA_UTOPIA_RANGE_HPP
#define UTOPIA_UTOPIA_RANGE_HPP


#include "utopia_Base.hpp"
#include <cassert>
#include <ostream>
#include <algorithm>


namespace utopia {
    class Range {
    private:
        const SizeType _begin, _end, _extent;

    public:
        Range(const SizeType begin, const SizeType to)
                : _begin(begin), _end(to), _extent(to - begin) { }


        explicit Range(const SizeType beginAndTo)
                : _begin(beginAndTo), _end(beginAndTo + 1), _extent(1)
        {
            assert(beginAndTo >= 0);
        }

        /** \addtogroup ranges
         *  @{
         */

        /*!
         * @return beginning of the range
         */
        inline SizeType begin() const
        {
            return _begin;
        }

        /*!
         * @return ending of the range. Hence, larger index contained in the range + 1
         */
        inline SizeType end() const
        {
            return _end;
        }

        /**
         * @return extent of the range => number of elements between 1st and the last element in the range. 
         */
        inline SizeType extent() const
        {
            return _extent;
        }

        /**
         * @brief      Checks if range is empty. 
         */
        inline bool empty() const {
            return _extent == 0;
        }

        /**
         * @brief      Checks if range is valid. 
         */
        inline bool valid() const {
            return _extent >= 0;
        }

        /**
         * @brief      Checks if given index is inside of the range. 
         */
        inline bool inside(const SizeType index) const
        {
            return index >= _begin && index < _end;
        }

        /**
         * @brief      Unites with other range.
         */
        inline Range unite(const Range &other) const {
            using std::max;
            using std::min;
            assert(are_contiguous(*this, other));
            return Range(min(begin(), other.begin()),
                         max(end(), other.end()));
        }

        /**
         * @brief      Finds intersection with other range. 
         */
        inline Range intersect(const Range &other) const {
            using std::max;
            using std::min;

            return Range(max(begin(), other.begin()),
                         min(end(), other.end()));
        }

        /*!
         * @return true of range1 and range2 || range2 and range1 are contiguos
         */
        inline friend bool are_contiguous(const Range &range1, const Range &range2)
        {
            return range1.begin() == range2.end() || range1.begin() == range2.end();
        }

        inline const Range operator + (const long shift) const {
            return Range(_begin + shift, _end + shift);
        }

        inline const Range operator - (const long shift) const {
            return Range(_begin - shift, _end - shift);
        }

        friend std::ostream & operator <<(std::ostream &os, const Range &range)
        {
            os << range.begin() << "\t";
            os << range.end() << "\t";
            os << "\n";
            return os;
        }

    /** @}*/

    };
}

#endif //UTOPIA_UTOPIA_RANGE_HPP

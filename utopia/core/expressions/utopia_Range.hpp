#ifndef UTOPIA_UTOPIA_RANGE_HPP
#define UTOPIA_UTOPIA_RANGE_HPP


#include "utopia_Base.hpp"
#include <cassert>
#include <ostream>
#include <algorithm>


namespace utopia {
    class Range {
    private:
        SizeType begin_, end_, extent_;

    public:
        Range() {}
        Range(const SizeType begin, const SizeType to)
                : begin_(begin), end_(to), extent_(to - begin) { }


        explicit Range(const SizeType beginAndTo)
                : begin_(beginAndTo), end_(beginAndTo + 1), extent_(1)
        {
            assert(beginAndTo >= 0);
        }

        inline void set(const SizeType &begin, const SizeType &end)
        {
            begin_ = begin;
            end_ = end;
            extent_ = end - begin;
        }

        /** \addtogroup ranges
         *  @{
         */

        /*!
         * @return beginning of the range
         */
        UTOPIA_INLINE_FUNCTION SizeType begin() const
        {
            return begin_;
        }

        /*!
         * @return ending of the range. Hence, larger index contained in the range + 1
         */
        UTOPIA_INLINE_FUNCTION SizeType end() const
        {
            return end_;
        }

        /**
         * @return extent of the range => number of elements between 1st and the last element in the range.
         */
        UTOPIA_INLINE_FUNCTION SizeType extent() const
        {
            return extent_;
        }

        /**
         * @brief      Checks if range is empty.
         */
        UTOPIA_INLINE_FUNCTION bool empty() const {
            return extent_ == 0;
        }

        /**
         * @brief      Checks if range is valid.
         */
        UTOPIA_INLINE_FUNCTION bool valid() const {
            return extent_ >= 0;
        }

        /**
         * @brief      Checks if given index is inside of the range.
         */
        UTOPIA_INLINE_FUNCTION bool inside(const SizeType index) const
        {
            return index >= begin_ && index < end_;
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
            return Range(begin_ + shift, end_ + shift);
        }

        inline const Range operator - (const long shift) const {
            return Range(begin_ - shift, end_ - shift);
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

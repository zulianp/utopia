#ifndef UTOPIA_UTOPIA_RANGE_HPP
#define UTOPIA_UTOPIA_RANGE_HPP

#include <algorithm>
#include <cassert>
#include <ostream>
#include "utopia_Algorithms.hpp"
#include "utopia_Base.hpp"

namespace utopia {

    class Range {
    private:
        SizeType begin_, end_, extent_;

    public:
        Range() = default;
        UTOPIA_INLINE_FUNCTION constexpr Range(const SizeType begin, const SizeType to)
            : begin_(begin), end_(to), extent_(to - begin) {}

        UTOPIA_INLINE_FUNCTION constexpr explicit Range(const SizeType beginAndTo)
            : begin_(beginAndTo), end_(beginAndTo + 1), extent_(1) {
            UTOPIA_DEVICE_ASSERT_CXX14(beginAndTo >= 0);
        }

        UTOPIA_INLINE_FUNCTION void set(const SizeType &begin, const SizeType &end) {
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
        UTOPIA_INLINE_FUNCTION constexpr SizeType begin() const { return begin_; }

        /*!
         * @return ending of the range. Hence, larger index contained in the range + 1
         */
        UTOPIA_INLINE_FUNCTION constexpr SizeType end() const { return end_; }

        /**
         * @return extent of the range => number of elements between 1st and the last element in the range.
         */
        UTOPIA_INLINE_FUNCTION constexpr SizeType extent() const { return extent_; }

        /**
         * @brief      Checks if range is empty.
         */
        UTOPIA_INLINE_FUNCTION constexpr bool empty() const { return extent_ == 0; }

        /**
         * @brief      Checks if range is valid.
         */
        UTOPIA_INLINE_FUNCTION constexpr bool valid() const { return extent_ >= 0; }

        /**
         * @brief      Checks if given index is inside of the range.
         */
        UTOPIA_INLINE_FUNCTION constexpr bool inside(const SizeType index) const {
            return index >= begin_ && index < end_;
        }

        /**
         * @brief      Unites with other range.
         */
        UTOPIA_INLINE_FUNCTION constexpr Range unite(const Range &other) const {
            UTOPIA_DEVICE_ASSERT_CXX14(are_contiguous(*this, other));
            return Range(device::min(begin(), other.begin()), device::max(end(), other.end()));
        }

        /**
         * @brief      Finds intersection with other range.
         */
        UTOPIA_INLINE_FUNCTION constexpr Range intersect(const Range &other) const {
            return Range(device::max(begin(), other.begin()), device::min(end(), other.end()));
        }

        /*!
         * @return true of range1 and range2 || range2 and range1 are contiguos
         */
        UTOPIA_INLINE_FUNCTION friend constexpr bool are_contiguous(const Range &range1, const Range &range2) {
            return range1.begin() == range2.end() || range1.begin() == range2.end();
        }

        UTOPIA_INLINE_FUNCTION const constexpr Range operator+(const long shift) const {
            return Range(begin_ + shift, end_ + shift);
        }

        UTOPIA_INLINE_FUNCTION const constexpr Range operator-(const long shift) const {
            return Range(begin_ - shift, end_ - shift);
        }

        UTOPIA_INLINE_FUNCTION constexpr bool operator==(const Range &other) const {
            return begin_ == other.begin_ && end_ == other.end_;
        }

        friend std::ostream &operator<<(std::ostream &os, const Range &range) {
            os << range.begin() << "\t";
            os << range.end() << "\t";
            os << "\n";
            return os;
        }

        /** @}*/
    };
}  // namespace utopia

#endif  // UTOPIA_UTOPIA_RANGE_HPP

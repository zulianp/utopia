#ifndef UTOPIA_ARRAY_VIEW_4_HPP
#define UTOPIA_ARRAY_VIEW_4_HPP

#include "utopia_ArrayView.hpp"

namespace utopia {

    template<typename T, Size_t N0_, Size_t N1_, Size_t N2_, Size_t N3_>
    class ArrayView<T, N0_, N1_, N2_, N3_> final {
    public:
        using SizeType = Size_t;
        static const SizeType Size = N0_ * N1_ * N2_ * N3_;

        UTOPIA_INLINE_FUNCTION constexpr SizeType size() const
        {
            return Size;
        }

        UTOPIA_INLINE_FUNCTION T &operator[](const SizeType &i)
        {
            UTOPIA_DEVICE_ASSERT(i < Size);

            return data_[i];
        }

        UTOPIA_INLINE_FUNCTION  T &operator()(
            const SizeType &i,
            const SizeType &j,
            const SizeType &k,
            const SizeType &l)
        {
            UTOPIA_DEVICE_ASSERT(i < N0_);
            UTOPIA_DEVICE_ASSERT(j < N1_);
            UTOPIA_DEVICE_ASSERT(k < N2_);
            UTOPIA_DEVICE_ASSERT(l < N3_);

            return data_[((i * N1_ + j) * N2_ + k) * N3_ + l];
        }

        UTOPIA_INLINE_FUNCTION constexpr const T &operator[](const SizeType &i) const
        {
            UTOPIA_DEVICE_ASSERT(i < Size);
            return data_[i];
        }

        UTOPIA_INLINE_FUNCTION constexpr const T &operator()(
            const SizeType &i,
            const SizeType &j,
            const SizeType &k,
            const SizeType &l) const
        {
            UTOPIA_DEVICE_ASSERT(i < N0_);
            UTOPIA_DEVICE_ASSERT(j < N1_);
            UTOPIA_DEVICE_ASSERT(k < N2_);
            UTOPIA_DEVICE_ASSERT(l < N3_);

            return (*this)[((i * N1_ + j) * N2_ + k) * N3_ + l];
        }

        UTOPIA_INLINE_FUNCTION constexpr SizeType extent(const SizeType &i) const
        {
            UTOPIA_DEVICE_ASSERT(i < 4);
            UTOPIA_DEVICE_ASSERT(i >= 0);

            switch(i) {
                case 0: return N0_;
                case 1: return N1_;
                case 2: return N2_;
                case 3: return N3_;
                default:
                {
                    return 0;
                }
            }
        }

        UTOPIA_INLINE_FUNCTION T* begin()
        {
            return data_.begin();
        }

        UTOPIA_INLINE_FUNCTION T* end()
        {
            return data_.end();
        }

        UTOPIA_INLINE_FUNCTION const T* begin() const
        {
            return data_.begin();
        }

        UTOPIA_INLINE_FUNCTION const T* end() const
        {
            return data_.end();
        }

        UTOPIA_INLINE_FUNCTION void copy(std::initializer_list<T> data)
        {
            SizeType i = 0;
            for(auto it = std::begin(data); it != std::end(data); ++it) {
                UTOPIA_DEVICE_ASSERT(i < Size);
                data_[i++] = *it;
            }
        }

        ArrayView(const ArrayView &other)
        {
            device::copy(other.data_, data_);
        }

        ArrayView &operator=(const ArrayView &other)
        {
            if(this == &other) return *this;
            device::copy(other.data_, data_);
            return *this;
        }

        ArrayView(ArrayView &&other)
        : data_(std::move(other.data_))
        {}

        ArrayView() {}

    private:
        ArrayView<T, Size> data_;
    };


}

#endif //UTOPIA_ARRAY_VIEW_4_HPP

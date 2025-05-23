#ifndef UTOPIA_ARRAY_2D_HPP
#define UTOPIA_ARRAY_2D_HPP

#include <array>
#include <string>

#include "utopia_Algorithms.hpp"
#include "utopia_Base.hpp"
#include "utopia_ForwardDeclarations.hpp"
#include "utopia_ViewForwardDeclarations.hpp"

namespace utopia {

    template <typename T, Size_t... Args>
    class ArrayView {};

    template <typename T>
    class ArrayView<T> {
    public:
        using SizeType = Size_t;

        virtual ~ArrayView() = default;

        UTOPIA_INLINE_FUNCTION constexpr SizeType size() const { return size_; }

        UTOPIA_INLINE_FUNCTION T &operator[](const SizeType &i) { return data_[i]; }

        UTOPIA_INLINE_FUNCTION constexpr const T &operator[](const SizeType &i) const { return data_[i]; }

        UTOPIA_INLINE_FUNCTION T *begin() { return data_; }

        UTOPIA_INLINE_FUNCTION T *end() { return data_ + size_; }

        UTOPIA_INLINE_FUNCTION const T *begin() const { return data_; }

        UTOPIA_INLINE_FUNCTION const T *end() const { return data_ + size_; }

        UTOPIA_INLINE_FUNCTION bool is_null() const { return data_ == nullptr; }

        UTOPIA_FUNCTION ArrayView() : data_(nullptr) {}
        UTOPIA_FUNCTION ArrayView(T *data, const SizeType &size) : data_(data), size_(size) {}

        template <typename OtherT, Size_t... SizeArgs>
        UTOPIA_INLINE_FUNCTION void copy(const ArrayView<OtherT, SizeArgs...> &other) {
            UTOPIA_DEVICE_ASSERT_CXX14(data_);
            UTOPIA_DEVICE_ASSERT_CXX14(size_ >= other.size());

            size_ = other.size();
            for (Size_t i = 0; i < size_; ++i) {
                (*this)[i] = other[i];
            }
        }

        template <Size_t... SizeArgs>
        UTOPIA_FUNCTION ArrayView(const ArrayView<T, SizeArgs...> &other) : data_(other.data_), size_(other.size()) {}

        UTOPIA_FUNCTION ArrayView(ArrayView &&other) : data_(std::move(other.data_)), size_(std::move(other.size_)) {}

        UTOPIA_FUNCTION ArrayView &operator=(ArrayView &&other) {
            data_ = std::move(other.data_);
            size_ = std::move(other.size_);
            return *this;
        }

        UTOPIA_FUNCTION void set(T *data, Size_t size) {
            data_ = data;
            size_ = size;
        }

        ArrayView &operator=(const ArrayView &other) = default;
        ArrayView(const ArrayView &other) = default;

        inline void set_data(T *data) { data_ = data; }
        UTOPIA_INLINE_FUNCTION bool empty() const { return size_ == 0; }

        template <typename OtherT>
        ArrayView<OtherT> reinterpret() {
            return ArrayView<OtherT>(reinterpret_cast<OtherT *>(data_), size_);
        }

    private:
        T *data_;
        Size_t size_{0};
    };

    template <typename T>
    class ArrayView<T, DYNAMIC_SIZE> final : public ArrayView<T> {};

    template <typename T, Size_t Size_>
    class ArrayView<T, Size_> final {
    public:
        using SizeType = Size_t;
        static const SizeType Size = Size_;

        UTOPIA_INLINE_FUNCTION constexpr SizeType size() const { return Size; }

        UTOPIA_INLINE_FUNCTION T &operator[](const SizeType &i) {
            UTOPIA_DEVICE_ASSERT_CXX14(i < size());

            return data_[i];
        }

        UTOPIA_INLINE_FUNCTION constexpr const T &operator[](const SizeType &i) const {
            UTOPIA_DEVICE_ASSERT_CXX14(i < size());

            return data_[i];
        }

        UTOPIA_INLINE_FUNCTION T *begin() { return data_; }

        UTOPIA_INLINE_FUNCTION T *end() { return data_ + Size; }

        UTOPIA_INLINE_FUNCTION constexpr const T *begin() const { return data_; }

        UTOPIA_INLINE_FUNCTION constexpr const T *end() const { return data_ + Size; }

        template <class ArrayViewOther>
        UTOPIA_INLINE_FUNCTION ArrayView &operator=(const ArrayViewOther &other) {
            // pragma GCCunroll(Size)
            for (Size_t i = 0; i < Size; ++i) {
                data_[i] = other[i];
            }

            return *this;
        }

        static UTOPIA_CONSTEXPR ArrayView make(const T &value) {
            ArrayView ret;

            // pragma GCCunroll(Size)
            for (Size_t i = 0; i < Size; ++i) {
                ret.data_[i] = value;
            }

            return ret;
        }

        ////public for aggregate intialization
        T data_[Size];
    };

    template <typename T, Size_t Rows_, Size_t Cols_>
    class ArrayView<T, Rows_, Cols_> final {
    public:
        using SizeType = Size_t;
        static const SizeType Rows = Rows_;
        static const SizeType Cols = Cols_;
        static const SizeType Size = Rows_ * Cols_;

        UTOPIA_INLINE_FUNCTION constexpr SizeType size() const { return Size; }

        UTOPIA_INLINE_FUNCTION T &operator[](const SizeType &i) {
            UTOPIA_DEVICE_ASSERT_CXX14(i < Size);

            return data_[i];
        }

        UTOPIA_INLINE_FUNCTION T &operator()(const SizeType &i, const SizeType &j) {
            UTOPIA_DEVICE_ASSERT_CXX14(i < Rows_);
            UTOPIA_DEVICE_ASSERT_CXX14(j < Cols_);

            return data_[i * Cols_ + j];
        }

        UTOPIA_INLINE_FUNCTION constexpr const T &operator[](const SizeType &i) const {
            UTOPIA_DEVICE_ASSERT_CXX14(i < Size);
            return data_[i];
        }

        UTOPIA_INLINE_FUNCTION constexpr const T &operator()(const SizeType &i, const SizeType &j) const {
            UTOPIA_DEVICE_ASSERT_CXX14(i < Rows_);
            UTOPIA_DEVICE_ASSERT_CXX14(j < Cols_);

            return data_[i * Cols_ + j];
        }

        UTOPIA_INLINE_FUNCTION constexpr SizeType extent(const SizeType &i) const {
            UTOPIA_DEVICE_ASSERT_CXX14(i < 2);
            UTOPIA_DEVICE_ASSERT_CXX14(i >= 0);

            return (i == 0) ? Rows_ : Cols_;
        }

        UTOPIA_INLINE_FUNCTION T *begin() { return data_.begin(); }

        UTOPIA_INLINE_FUNCTION T *end() { return data_.end(); }

        UTOPIA_INLINE_FUNCTION const T *begin() const { return data_.begin(); }

        UTOPIA_INLINE_FUNCTION const T *end() const { return data_.end(); }

        UTOPIA_INLINE_FUNCTION void copy(std::initializer_list<T> data) {
            SizeType i = 0;
            for (auto it = std::begin(data); it != std::end(data); ++it) {
                UTOPIA_DEVICE_ASSERT_CXX14(i < Size);
                data_[i++] = *it;
            }
        }

        ArrayView(const ArrayView &other) { device::copy(other.data_, data_); }

        ArrayView &operator=(const ArrayView &other) {
            if (this == &other) return *this;
            device::copy(other.data_, data_);
            return *this;
        }

        ArrayView(ArrayView &&other) : data_(std::move(other.data_)) {}

        ArrayView() = default;

    private:
        ArrayView<T, Size> data_;
    };

    template <typename T>
    class ArrayView<T, DYNAMIC_SIZE, DYNAMIC_SIZE> final {
    public:
        using SizeType = Size_t;

        UTOPIA_INLINE_FUNCTION constexpr SizeType size() const { return rows_ * cols_; }

        UTOPIA_INLINE_FUNCTION T &operator[](const SizeType &i) {
            UTOPIA_DEVICE_ASSERT_CXX14(i < size());

            return data_[i];
        }

        UTOPIA_INLINE_FUNCTION T &operator()(const SizeType &i, const SizeType &j) {
            UTOPIA_DEVICE_ASSERT_CXX14(i < extent(0));
            UTOPIA_DEVICE_ASSERT_CXX14(j < extent(1));

            return data_[i * cols_ + j];
        }

        UTOPIA_INLINE_FUNCTION constexpr const T &operator[](const SizeType &i) const {
            UTOPIA_DEVICE_ASSERT_CXX14(i < size());

            return data_[i];
        }

        UTOPIA_INLINE_FUNCTION constexpr const T &operator()(const SizeType &i, const SizeType &j) const {
            UTOPIA_DEVICE_ASSERT_CXX14(i < extent(0));
            UTOPIA_DEVICE_ASSERT_CXX14(j < extent(1));

            return data_[i * cols_ + j];
        }

        UTOPIA_INLINE_FUNCTION constexpr SizeType extent(const SizeType &i) const {
            UTOPIA_DEVICE_ASSERT_CXX14(i < 2);
            UTOPIA_DEVICE_ASSERT_CXX14(i >= 0);

            return (i == 0) ? rows_ : cols_;
        }

        UTOPIA_INLINE_FUNCTION T *begin() { return data_.begin(); }

        UTOPIA_INLINE_FUNCTION T *end() { return data_.end(); }

        UTOPIA_INLINE_FUNCTION const T *begin() const { return data_.begin(); }

        UTOPIA_INLINE_FUNCTION const T *end() const { return data_.end(); }

        UTOPIA_FUNCTION ArrayView(T *data, const SizeType &rows, const SizeType &cols)
            : data_(data, rows * cols), rows_(rows), cols_(cols) {}

        UTOPIA_FUNCTION ArrayView() : data_() {}

        UTOPIA_INLINE_FUNCTION bool is_null() const { return data_.is_null(); }

        inline void set_data(T *data) { data_.set_data(data); }
        inline void set_size(const SizeType rows, const SizeType cols) {
            rows_ = rows;
            cols_ = cols;
        }

    private:
        ArrayView<T> data_;
        SizeType rows_{0}, cols_{0};
    };

    template <typename T, Size_t... Args>
    void disp(const ArrayView<T, Args...> &view, std::ostream &os = std::cout) {
        for (const auto &v : view) {
            os << v << " ";
        }

        os << "\n";
    }

    template <typename T, size_t Size>
    inline void convert(const ArrayView<T, Size> &in, std::array<T, Size> &out) {
        std::copy(&in[0], &in[0] + Size, &out[0]);
    }

}  // namespace utopia

#endif  // UTOPIA_ARRAY_2D_HPP

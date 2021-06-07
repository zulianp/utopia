#ifndef UTOPIA_BLAS_ARRAY_HPP
#define UTOPIA_BLAS_ARRAY_HPP

#include <vector>
#include "utopia_Array.hpp"

namespace utopia {

    template <typename EntryType_>
    class BlasArray final : public Array<EntryType_, std::size_t> {
    public:
        using EntryType = EntryType_;
        using SizeType = std::size_t;
        using Iterator = typename std::vector<EntryType>::iterator;
        using ConstIterator = typename std::vector<EntryType>::const_iterator;

        ~BlasArray() {}

        // locks
        inline void read_lock() override {}
        inline void write_lock(WriteMode) override {}

        inline void read_unlock() override {}
        inline void write_unlock(WriteMode) override {}

        // basic mutators
        inline void set(const SizeType &i, const EntryType &value) override {
            assert(i < size());
            entries_[i] = value;
        }

        inline EntryType get(const SizeType &i) override {
            assert(i < size());
            return entries_[i];
        }

        inline EntryType &operator[](const SizeType &i) {
            assert(i < size());
            return entries_[i];
        }

        inline const EntryType &operator[](const SizeType &i) const {
            assert(i < size());
            return entries_[i];
        }

        // print function
        inline void describe() const override {
            for (auto e : entries_) {
                utopia::out() << e << " ";
            }

            utopia::out() << std::endl;
        }

        // utility functions
        inline bool empty() const override { return entries_.empty(); }

        inline void clear() override { entries_.clear(); }

        inline SizeType size() const override { return entries_.size(); }

        void resize(const SizeType size) { entries_.resize(size); }

        inline Iterator begin() { return entries_.begin(); }

        inline ConstIterator begin() const { return entries_.begin(); }

        inline Iterator end() { return entries_.end(); }

        inline ConstIterator end() const { return entries_.end(); }

    private:
        std::vector<EntryType> entries_;
    };

}  // namespace utopia

#endif  // UTOPIA_BLAS_ARRAY_HPP

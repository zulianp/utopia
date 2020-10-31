#ifndef UTOPIA_BLAS_INDEX_SET_HPP
#define UTOPIA_BLAS_INDEX_SET_HPP

#include <utility>
#include <vector>
#include "utopia_IndexSet.hpp"

namespace utopia {

    class BlasIndexSet : public IndexSet<int, std::size_t> {
    public:
        using Index = int;
        using SizeType = std::size_t;

        template <class... Args>
        BlasIndexSet(Args &&... args) : index_(std::forward<Args>(args)...) {}

        ~BlasIndexSet() override = default;

        // locks
        inline void read_lock() override {}
        inline void write_lock(WriteMode) override {}

        inline void read_unlock() override {}
        inline void write_unlock(WriteMode) override {}

        // basic mutators
        void set(const SizeType &i, const Index &value) override {
            assert(i < size());
            index_[i] = value;
        }

        Index get(const SizeType &i) const override {
            assert(i < size());
            return index_[i];
        }

        Index &operator[](const SizeType &i) {
            assert(i < size());
            return index_[i];
        }

        const Index &operator[](const SizeType &i) const {
            assert(i < size());
            return index_[i];
        }

        // print function
        inline void describe() const override {
            for (auto i : index_) {
                utopia::out() << i << " ";
            }

            utopia::out() << std::endl;
        }

        // utility functions
        inline bool empty() const override { return index_.empty(); }

        inline void clear() override { index_.clear(); }

        inline void resize(const SizeType n) { index_.resize(n); }

        inline SizeType size() const override { return index_.size(); }

    private:
        std::vector<Index> index_;
    };

}  // namespace utopia

#endif  // UTOPIA_BLAS_INDEX_SET_HPP

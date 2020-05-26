#ifndef UTOPIA_ELEM_HPP
#define UTOPIA_ELEM_HPP

#include <string>

namespace utopia {

    class Elem {
    public:
        Elem() = default;
        virtual ~Elem() = default;

        const std::size_t &idx() const { return idx_; }

        void idx(const std::size_t &idx) { idx_ = idx; }

    private:
        std::size_t idx_{static_cast<std::size_t>(-1)};
    };
}  // namespace utopia

#endif  // UTOPIA_ELEM_HPP

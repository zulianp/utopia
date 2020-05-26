#ifndef UTOPIA_RENAME_HPP
#define UTOPIA_RENAME_HPP

#include "utopia_Traits.hpp"

namespace utopia {
    template <class Tensor, int Backend = Traits<Tensor>::Backend>
    class Rename {
    public:
        inline static void apply(const std::string &, const Tensor &) {}
    };

    template <class Tensor>
    void rename(const std::string &name, Tensor &t) {
        Rename<Tensor>::apply(name, t);
    }
}  // namespace utopia

#endif  // UTOPIA_RENAME_HPP

#ifndef UTOPIA_PETSC_EVAL_RENAME_HPP
#define UTOPIA_PETSC_EVAL_RENAME_HPP

#include "utopia_Rename.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    template <class Tensor>
    class Rename<Tensor, PETSC> {
    public:
        inline static void apply(const std::string &name, Tensor &t) { t.rename(name); }
    };
}  // namespace utopia

#endif  // UTOPIA_PETSC_EVAL_RENAME_HPP

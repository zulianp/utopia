#ifndef UTOPIA_GLOBAL_INDEX_HPP
#define UTOPIA_GLOBAL_INDEX_HPP

#include "utopia_ViewForwardDeclarations.hpp"

namespace utopia {

    template <class FunctionSpace, int BlockSize = utopia::DYNAMIC_SIZE>
    class GlobalIndex {
    public:
        GlobalIndex(int n_var) : n_var_(n_var) {}

    private:
        int n_var_;
    };

}  // namespace utopia

#endif  // UTOPIA_GLOBAL_INDEX_HPP
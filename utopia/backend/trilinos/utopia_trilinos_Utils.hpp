#ifndef UTOPIA_TRILINOS_UTILS_HPP
#define UTOPIA_TRILINOS_UTILS_HPP

#include "utopia_Communicator.hpp"

namespace utopia {

    template <typename SizeType>
    inline SizeType decompose(const Communicator &comm, const SizeType n_global) {
        const SizeType rank = comm.rank();
        const SizeType size = comm.size();

        const SizeType n = n_global / size;
        const SizeType reminder = n_global % size;
        const SizeType n_local = n + static_cast<SizeType>(rank < reminder);
        return n_local;
    }
}  // namespace utopia

#endif  // UTOPIA_TRILINOS_UTILS_HPP

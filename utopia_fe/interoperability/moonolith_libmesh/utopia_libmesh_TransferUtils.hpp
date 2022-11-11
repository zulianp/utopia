#ifndef UTOPIA_TRANSFER_UTILS_HPP
#define UTOPIA_TRANSFER_UTILS_HPP

#include "utopia_fe_base.hpp"

namespace utopia {

    void tensorize(const USparseMatrix &T_x, const SizeType n_var, USparseMatrix &T);
    void tensorize(const SizeType n_var, UVector &t);

}  // namespace utopia

#endif  // UTOPIA_TRANSFER_UTILS_HPP

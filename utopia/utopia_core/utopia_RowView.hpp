#ifndef UTOPIA_ROW_VIEW_HPP
#define UTOPIA_ROW_VIEW_HPP

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"

namespace utopia {
    template <class Tensor,
              int Order = Tensor::Order,
              int FILL_TYPE = Traits<Tensor>::FILL_TYPE,
              int BackendType = Traits<Tensor>::Backend>
    class RowView {};
}  // namespace utopia

#endif  // UTOPIA_ROW_VIEW_HPP

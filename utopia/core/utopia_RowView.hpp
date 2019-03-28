#ifndef UTOPIA_ROW_VIEW_HPP
#define UTOPIA_ROW_VIEW_HPP

#include "utopia_Base.hpp"

namespace utopia {
    template<class Tensor, int Order = Tensor::Order, int FILL_TYPE = Tensor::FILL_TYPE, int BackendType = Traits<Tensor>::Backend>
    class RowView {};
}

#endif //UTOPIA_ROW_VIEW_HPP


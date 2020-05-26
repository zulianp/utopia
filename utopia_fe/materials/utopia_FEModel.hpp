#ifndef UTOPIA_FE_MODEL_HPP
#define UTOPIA_FE_MODEL_HPP

#include "utopia_Model.hpp"

namespace utopia {

    template <class FunctionSpace, class Matrix, class Vector>
    class FEModel : public Model<Matrix, Vector> {
    public:
        virtual ~FEModel() {}

        virtual FunctionSpace &space() = 0;
        virtual const FunctionSpace &space() const = 0;
    };
}  // namespace utopia

#endif  // UTOPIA_FE_MODEL_HPP

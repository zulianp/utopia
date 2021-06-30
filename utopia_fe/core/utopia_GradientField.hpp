#ifndef UTOPIA_GRADIENT_FIELD_HPP
#define UTOPIA_GRADIENT_FIELD_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_Instance.hpp"

#include "utopia_Field.hpp"

#include <memory>
#include <string>

namespace utopia {

    template <class FunctionSpace>
    class GradientField : public Field<FunctionSpace> {
    public:
        using Vector = typename Traits<FunctionSpace>::Vector;

        bool init(Field<FunctionSpace> &field) { return false; }

        void normalize() {}
    };

}  // namespace utopia

#endif  // UTOPIA_GRADIENT_FIELD_HPP

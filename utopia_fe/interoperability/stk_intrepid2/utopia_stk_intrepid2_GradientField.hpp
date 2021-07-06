#ifndef UTOPIA_STK_INTREPID2_GRADIENTFIELD_HPP
#define UTOPIA_STK_INTREPID2_GRADIENTFIELD_HPP

#include "utopia_Input.hpp"
#include "utopia_Traits.hpp"

#include "utopia_Instance.hpp"

#include "utopia_Field.hpp"
#include "utopia_GradientField.hpp"

#include "utopia_stk_FunctionSpace.hpp"

#include <memory>
#include <string>

namespace utopia {

    template <>
    class GradientField<stk::FunctionSpace> : public Field<stk::FunctionSpace> {
    public:
        using Vector = typename Traits<stk::FunctionSpace>::Vector;

        bool init_and_normalize(Field<stk::FunctionSpace> &field);

        void normalize();
    };

}  // namespace utopia

#endif  // UTOPIA_STK_INTREPID2_GRADIENTFIELD_HPP
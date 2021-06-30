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

        bool init(Field<FunctionSpace> &field) {
            this->set_space(field.space());
            this->set_name("grad_" + field.name());

            auto dim = field.space()->mesh().spatial_dimension();
            this->set_offset(field.offset());
            this->set_tensor_size(field.tensor_size() * dim);

            auto l = layout(field.data());
            auto lg = layout(l.comm(), l.local_size(), l.size());
            auto data = std::make_shared<Vector>(lg);

            // utopia::grad(field, *this);

            this->set_data(data);

            return false;
        }

        void normalize() {}
    };

}  // namespace utopia

#endif  // UTOPIA_GRADIENT_FIELD_HPP

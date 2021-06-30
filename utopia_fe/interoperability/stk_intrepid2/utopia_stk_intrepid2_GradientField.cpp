#include "utopia_stk_intrepid2_GradientField.hpp"

#include "utopia_stk_intrepid2.hpp"

// Utopia/intrepid2 includes
#include "utopia_intrepid2.hpp"

namespace utopia {

    bool GradientField<stk::FunctionSpace>::init(Field<stk::FunctionSpace> &field) {
        using Scalar_t = Traits<stk::FunctionSpace>::Scalar;
        using Intrepid2Field_t = utopia::intrepid2::Field<Scalar_t>;
        using Intrepid2Gradient_t = utopia::intrepid2::Gradient<Scalar_t>;
        using Intrepid2FE_t = utopia::intrepid2::FE<Scalar_t>;

        this->set_space(field.space());
        this->set_name("grad_" + field.name());

        auto dim = field.space()->mesh().spatial_dimension();
        this->set_offset(field.offset());
        this->set_tensor_size(field.tensor_size() * dim);

        auto l = layout(field.data());
        auto lg = layout(l.comm(), l.local_size(), l.size());
        auto data = std::make_shared<Vector>(lg);

        auto fe = std::make_shared<Intrepid2FE_t>();

        // One gradient per element
        int quadrature_order = 2;
        create_fe(*this->space(), *fe, quadrature_order);

        Intrepid2Field_t intrepid_field(fe);
        Intrepid2Gradient_t intrepid_gradient(fe);

        convert_field(field, intrepid_field);
        intrepid_gradient.init(intrepid_field);

        // Compute mass contributions and projections

        this->set_data(data);

        return false;
    }

    void GradientField<stk::FunctionSpace>::normalize() {}

}  // namespace utopia

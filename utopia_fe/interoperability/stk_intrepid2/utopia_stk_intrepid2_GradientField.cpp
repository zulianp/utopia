#include "utopia_stk_intrepid2_GradientField.hpp"

#include "utopia_kokkos_Commons.hpp"
#include "utopia_kokkos_Gradient.hpp"

// Utopia/intrepid2 includes
#include "utopia_intrepid2.hpp"

// Utopia/Stk/intrepid2 includes
#include "utopia_stk_intrepid2.hpp"
#include "utopia_stk_intrepid2_L2Projection.hpp"
#include "utopia_stk_intrepid2_OmniAssembler.hpp"

namespace utopia {

    bool GradientField<stk::FunctionSpace>::init_and_normalize(Field<stk::FunctionSpace> &field) {
        using Scalar_t = Traits<stk::FunctionSpace>::Scalar;

        using Intrepid2FE_t = utopia::intrepid2::FE<Scalar_t>;
        using Intrepid2Gradient_t = utopia::kokkos::Gradient<Intrepid2FE_t>;
        using Intrepid2Field_t = utopia::kokkos::Field<Intrepid2FE_t>;

        this->set_space(field.space());
        this->set_name("grad_" + field.name());

        auto dim = field.space()->mesh().spatial_dimension();
        this->set_offset(field.offset());
        this->set_tensor_size(field.tensor_size() * dim);

        auto l = layout(field.data());
        auto lg = layout(l.comm(), l.local_size() * dim, l.size() * dim);
        auto data = std::make_shared<Vector>(lg);

        auto fe = std::make_shared<Intrepid2FE_t>();

        int quadrature_order = 0;
        create_fe(*this->space(), *fe, quadrature_order);

        Intrepid2Field_t intrepid_field(fe);
        Intrepid2Gradient_t intrepid_gradient(fe);

        // global_to_local(*field.space(), field.data(), intrepid_field.data());
        convert_field(field, intrepid_field);

        intrepid_gradient.init(intrepid_field);

        // Compute mass contributions and projections
        stk::L2Projection assembler(fe);
        assembler.set_field(make_ref(intrepid_gradient));
        assembler.init();
        assembler.ensure_vector_accumulator();
        if (!assembler.assemble_vector()) {
            return false;
        }

        data->set(0.0);
        assert(assembler.vector_accumulator());
        auto el_vec = assembler.vector_data();

        local_to_global(*field.space(), el_vec, ADD_MODE, *data, dim);

        this->set_data(data);

        normalize();
        return true;
    }

    void GradientField<stk::FunctionSpace>::normalize() {
        using Scalar_t = Traits<stk::FunctionSpace>::Scalar;

        auto view = local_view_device(this->data());

        auto dim = this->tensor_size();

        auto rd_complete = local_range_device(this->data());
        RangeDevice<Vector> rd(rd_complete.begin(), rd_complete.begin() + rd_complete.extent() / dim);

        parallel_for(
            rd, UTOPIA_LAMBDA(const int i) {
                Scalar_t norm_v = 0.;
                for (int d = 0; d < dim; ++d) {
                    auto x = view.get(i * dim + d);
                    norm_v += x * x;
                }

                if (norm_v > 0) {
                    // assert(norm_v > 0);
                    norm_v = device::sqrt(norm_v);

                    for (int d = 0; d < dim; ++d) {
                        auto x = view.get(i * dim + d);
                        view.set(i * dim + d, x / norm_v);
                    }
                }
            });
    }

}  // namespace utopia

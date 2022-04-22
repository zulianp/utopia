#ifndef UTOPIA_INTREPI2_FLOW_POST_PROCESS_APP_HPP
#define UTOPIA_INTREPI2_FLOW_POST_PROCESS_APP_HPP

// Utopia includes
#include "utopia_CreateFE.hpp"
#include "utopia_Field.hpp"

// Utopia/Kokkos includes
#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_ForcingFunction.hpp"
#include "utopia_kokkos_Gradient.hpp"
#include "utopia_kokkos_LaplaceOperator.hpp"

// Utopia/intrepid2 includes
#include "utopia_intrepid2.hpp"

namespace utopia {
    namespace intrepid2 {

        template <class FunctionSpace>
        class FlowPostProcessApp : public Configurable {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using Vector_t = typename Traits<FunctionSpace>::Vector;
            using Field_t = utopia::Field<FunctionSpace>;
            using IO_t = utopia::IO<FunctionSpace>;

            using Intrepid2FE_t = utopia::intrepid2::FE<Scalar_t>;
            using Intrepid2Field_t = utopia::kokkos::Field<Intrepid2FE_t>;
            using Intrepid2QPField_t = utopia::kokkos::QPField<Intrepid2FE_t>;
            using Intrepid2Gradient_t = utopia::kokkos::Gradient<Intrepid2FE_t>;

            using Flow_t = utopia::kokkos::LaplaceOperator<Intrepid2FE_t>;
            using ForcingFunction_t = utopia::kokkos::ForcingFunction<Intrepid2FE_t>;

            void read(Input &in) override {
                valid_ = true;

                if (!Options()
                         .add_option("verbose", verbose, "Verbose output.")
                         .add_option("quadrature_order", quadrature_order, "Specify the quadrature order.")
                         // .add_option("output_path", output_path, "Path to output file.")
                         .parse(in)) {
                    valid_ = false;
                    return;
                }

                in.get("space", [this](Input &in) {
                    space.read_with_state(in, field);

                    if (verbose) {
                        const Scalar_t norm_field = norm2(field.data());
                        std::cout << "norm_field: " << norm_field << std::endl;
                    }
                });

                if (space.empty()) {
                    valid_ = false;
                    return;
                }

                // Intrepid 2 stuff
                fe = std::make_shared<Intrepid2FE_t>();
                create_fe(space, *fe, quadrature_order);

                material = std::make_shared<Flow_t>(fe);
                in.get("material", *material);
            }

            bool valid() const { return valid_; }

            void run() {
                if (!valid()) return;

                if (verbose) {
                    int dim = space.mesh().spatial_dimension();
                    utopia::out() << "Dim:" << dim << " \n";
                }

                ////////////////////////////////////////////////////

                // if (export_strains) {
                //     Field_t avg_principal_strains("principal_strains", make_ref(space),
                //     std::make_shared<Vector_t>());

                //     {
                //         // Intrepid 2 stuff
                //         Intrepid2Field_t i2_displacement(fe);
                //         Intrepid2Gradient_t i2_gradient(fe);
                //         Intrepid2Strain_t i2_intrepid_strain(fe);
                //         Intrepid2QPField_t i2_principal_strains(fe);
                //         Intrepid2Field_t i2_avg_principal_strains(fe);

                //         convert_field(field, i2_displacement);
                //         i2_gradient.init(i2_displacement);

                //         if (!material->is_linear()) {
                //             i2_gradient.add_identity();
                //             // auto intrepid_det = det(i2_gradient);
                //         } else {
                //             i2_intrepid_strain.init_linearized(i2_displacement);
                //         }

                //         i2_intrepid_strain.eig(i2_principal_strains);
                //         i2_principal_strains.avg(i2_avg_principal_strains);
                //         i2_avg_principal_strains.set_elem_type(ELEMENT_TYPE);

                //         convert_field(i2_avg_principal_strains, avg_principal_strains);
                //     }

                //     {
                //         space.displace(field.data());  // FIXME
                //         IO_t io(space);
                //         io.set_output_path(output_path);
                //         io.write(avg_principal_strains);
                //     }
                // }

                // ////////////////////////////////////////////////////

                // else if (material) {
                //     Field_t principal_stresses("principal_stresses", make_ref(space), std::make_shared<Vector_t>());

                //     {
                //         // Intrepid 2 stuff
                //         Intrepid2Field_t i2_displacement(fe);
                //         Intrepid2Field_t i2_principal_stresses(fe);

                //         convert_field(field, i2_displacement);

                //         material->principal_stresses(i2_displacement.data(), i2_principal_stresses.data());
                //         i2_principal_stresses.set_elem_type(ELEMENT_TYPE);

                //         convert_field(i2_principal_stresses, principal_stresses);
                //     }

                //     {
                //         space.displace(field.data());  // FIXME
                //         IO_t io(space);
                //         io.set_output_path(output_path);
                //         io.write(principal_stresses);
                //     }
                // }
            }

            FunctionSpace space;
            Field_t field;

            bool valid_{false};
            bool verbose{false};
            int quadrature_order{0};

            std::shared_ptr<Intrepid2FE_t> fe;
            std::shared_ptr<Flow_t> material;
            std::vector<std::shared_ptr<ForcingFunction_t>> forcing_functions;
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPI2_FLOW_POST_PROCESS_APP_HPP
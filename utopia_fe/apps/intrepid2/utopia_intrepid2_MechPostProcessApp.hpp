#ifndef UTOPIA_INTREPID2_MECH_POST_PROCESS_APP_HPP
#define UTOPIA_INTREPID2_MECH_POST_PROCESS_APP_HPP

// Utopia includes
#include "utopia_CreateFE.hpp"
#include "utopia_Field.hpp"

// Utopia/Kokkos includes
#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_Gradient.hpp"
#include "utopia_kokkos_Strain.hpp"

// Utopia/intrepid2 includes
#include "utopia_intrepid2.hpp"

namespace utopia {
    namespace intrepid2 {

        template <class FunctionSpace>
        class MechPostProcessApp : public Configurable {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using Vector_t = typename Traits<FunctionSpace>::Vector;
            using Field_t = utopia::Field<FunctionSpace>;
            using IO_t = utopia::IO<FunctionSpace>;

            using Intrepid2FE_t = utopia::intrepid2::FE<Scalar_t>;
            using Intrepid2Field_t = utopia::kokkos::Field<Intrepid2FE_t>;
            using Intrepid2QPField_t = utopia::kokkos::QPField<Intrepid2FE_t>;
            using Intrepid2Gradient_t = utopia::kokkos::Gradient<Intrepid2FE_t>;
            using Intrepid2Strain_t = utopia::kokkos::Strain<Intrepid2FE_t>;

            // Dummy
            class Material {
            public:
                bool is_linear() const { return true; }
            };

            void read(Input &in) override {
                valid_ = true;

                if (!Options()
                         .add_option("verbose", verbose, "Verbose output.")
                         .add_option("quadrature_order", quadrature_order, "Specify the quadrature order.")
                         .add_option("output_path", output_path, "Path to output file.")
                         .parse(in)) {
                    valid_ = false;
                    return;
                }

                in.get("space", [this](Input &in) {
                    space.read_with_state(in, displacement_field);

                    if (verbose) {
                        const Scalar_t norm_field = norm2(displacement_field.data());
                        std::cout << "norm_field: " << norm_field << std::endl;
                    }
                });

                if (space.empty()) {
                    valid_ = false;
                    return;
                }

                material = std::make_shared<Material>();
            }

            bool valid() const { return valid_; }

            void run() {
                if (!valid()) return;

                if (verbose) {
                    int dim = space.mesh().spatial_dimension();
                    utopia::out() << "Dim:" << dim << " \n";
                }

                Field_t avg_principal_strains("principal_strains", make_ref(space), std::make_shared<Vector_t>());

                ////////////////////////////////////////////////////

                {
                    // Intrepid 2 stuff

                    auto fe = std::make_shared<Intrepid2FE_t>();

                    // One gradient per element
                    create_fe(space, *fe, quadrature_order);

                    Intrepid2Field_t i2_displacement(fe);
                    Intrepid2Gradient_t i2_gradient(fe);
                    Intrepid2Strain_t i2_intrepid_strain(fe);
                    Intrepid2QPField_t i2_principal_strains(fe);
                    Intrepid2Field_t i2_avg_principal_strains(fe);

                    convert_field(displacement_field, i2_displacement);
                    i2_gradient.init(i2_displacement);

                    if (!material->is_linear()) {
                        i2_gradient.add_identity();
                        // auto intrepid_det = det(i2_gradient);
                    } else {
                        i2_intrepid_strain.init_linearized(i2_displacement);
                    }

                    i2_intrepid_strain.eig(i2_principal_strains);
                    i2_principal_strains.avg(i2_avg_principal_strains);
                    i2_avg_principal_strains.set_elem_type(ELEMENT_TYPE);

                    convert_field(i2_avg_principal_strains, avg_principal_strains);
                }

                ////////////////////////////////////////////////////

                IO_t io(space);
                io.set_output_path(output_path);
                io.write(avg_principal_strains);
            }

            FunctionSpace space;
            Field_t displacement_field;

            bool valid_{false};
            bool verbose{false};
            int quadrature_order{0};

            Path output_path{"out.e"};

            std::shared_ptr<Material> material;
        };

    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_MECH_POST_PROCESS_APP_HPP
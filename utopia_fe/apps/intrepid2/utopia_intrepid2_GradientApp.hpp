#ifndef UTOPIA_INTREPID2_GRADIENT_APP_HPP
#define UTOPIA_INTREPID2_GRADIENT_APP_HPP

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
        class GradientApp : public Configurable {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using Vector_t = typename Traits<FunctionSpace>::Vector;
            using Field_t = utopia::Field<FunctionSpace>;
            using IO_t = utopia::IO<FunctionSpace>;

            using Intrepid2FE_t = utopia::intrepid2::FE<Scalar_t>;
            using Intrepid2Field_t = utopia::kokkos::Field<Intrepid2FE_t>;
            using Intrepid2Gradient_t = utopia::kokkos::Gradient<Intrepid2FE_t>;
            using Intrepid2Strain_t = utopia::kokkos::Strain<Intrepid2FE_t>;

            void read(Input &in) override {
                valid_ = true;

                if (!Options()
                         .add_option("verbose", verbose_, "Verbose output.")
                         .add_option("compute_deformation_gradient",
                                     compute_deformation_gradient_,
                                     "Compute and analyse the deformation gradient.")
                         .add_option("quadrature_order", quadrature_order, "Specify the quadrature order.")
                         .add_option("compute_strain", compute_strain_, "Compute and analyse the strain tensor")
                         .parse(in)) {
                    valid_ = false;
                    return;
                }

                in.get("space", [this](Input &in) {
                    space_.read_with_state(in, field_);

                    if (verbose_) {
                        const Scalar_t norm_field = norm2(field_.data());
                        std::cout << "norm_field: " << norm_field << std::endl;
                    }
                });

                if (space_.empty()) {
                    valid_ = false;
                    return;
                }
            }

            bool valid() const { return valid_; }

            void run() {
                if (!valid()) return;

                auto fe = std::make_shared<Intrepid2FE_t>();

                // One gradient per element
                create_fe(space_, *fe, quadrature_order);

                Intrepid2Field_t intrepid_field(fe);
                Intrepid2Gradient_t intrepid_gradient(fe);

                convert_field(field_, intrepid_field);

                intrepid_gradient.init(intrepid_field);

                if (verbose_) {
                    utopia::out() << "Gradient:\n";
                    utopia::out() << "-------------------------------\n";
                    intrepid_gradient.describe(utopia::out().stream());
                    utopia::out() << "-------------------------------\n";
                }

                Intrepid2Field_t avg_gradient(fe);
                intrepid_gradient.avg(avg_gradient);
                avg_gradient.set_elem_type(ELEMENT_TYPE);

                Field_t avg_gradient_global("avggrad", make_ref(space_), std::make_shared<Vector_t>());
                convert_field(avg_gradient, avg_gradient_global);

                utopia::out() << avg_gradient_global.data().size() << "\n";
                utopia::out() << avg_gradient_global.tensor_size() << "\n";

                {
                    IO_t io(space_);
                    io.set_output_path(path_);
                    io.write(avg_gradient_global);
                }

                if (verbose_) {
                    utopia::out() << "avg(Gradient):\n";
                    utopia::out() << "-------------------------------\n";
                    avg_gradient.describe(utopia::out().stream());
                    utopia::out() << "-------------------------------\n";
                }

                if (compute_deformation_gradient_) {
                    intrepid_gradient.add_identity();

                    auto intrepid_det = det(intrepid_gradient);

                    if (verbose_) {
                        utopia::out() << "Deformation gradient:\n";
                        utopia::out() << "-------------------------------\n";
                        intrepid_gradient.describe(utopia::out().stream());
                        utopia::out() << "-------------------------------\n";
                    }

                    utopia::out() << "Deformation gradient determinant:\n";
                    intrepid_det.describe(utopia::out().stream());
                }

                if (compute_strain_) {
                    Intrepid2Strain_t intrepid_strain(fe);
                    intrepid_strain.init_linearized(intrepid_field);

                    Intrepid2Field_t avg_strain(fe);
                    intrepid_strain.avg(avg_strain);
                    avg_strain.set_elem_type(ELEMENT_TYPE);

                    Field_t avg_strain_global("avgstrain", make_ref(space_), std::make_shared<Vector_t>());
                    convert_field(avg_strain, avg_strain_global);

                    {
                        IO_t io(space_);
                        io.set_output_path("strain.e");
                        io.write(avg_strain_global);
                    }

                    if (verbose_) {
                        utopia::out() << "Strain:\n";
                        utopia::out() << "-------------------------------\n";
                        intrepid_strain.describe(utopia::out().stream());
                        utopia::out() << "-------------------------------\n";
                    }

                    auto intrepid_det = det(intrepid_strain);
                    utopia::out() << "Strain determinant:\n";
                    intrepid_det.describe(utopia::out().stream());
                }
            }

        private:
            FunctionSpace space_;
            Field_t field_;

            bool valid_{false};
            bool verbose_{false};
            bool compute_deformation_gradient_{false};
            bool compute_strain_{false};
            int quadrature_order{0};

            Path path_{"out.e"};
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_GRADIENT_APP_HPP

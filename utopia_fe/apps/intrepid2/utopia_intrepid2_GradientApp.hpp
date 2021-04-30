#ifndef UTOPIA_INTREPID2_GRADIENT_APP_HPP
#define UTOPIA_INTREPID2_GRADIENT_APP_HPP

// Utopia includes
#include "utopia_CreateFE.hpp"
#include "utopia_Field.hpp"

// Utopia/intrepid2 includes
#include "utopia_intrepid2.hpp"

namespace utopia {
    namespace intrepid2 {

        template <class FunctionSpace>
        class GradientApp : public Configurable {
        public:
            using Scalar_t = typename Traits<FunctionSpace>::Scalar;
            using Field_t = utopia::Field<FunctionSpace>;

            using Intrepid2Field_t = utopia::intrepid2::Field<Scalar_t>;
            using Intrepid2Gradient_t = utopia::intrepid2::Gradient<Scalar_t>;
            using Intrepid2FE_t = utopia::intrepid2::FE<Scalar_t>;

            void read(Input &in) override {
                valid_ = true;

                if (!Options()
                         .add_option("verbose", verbose_, "Verbose output.")
                         .add_option("compute_deformation_gradient",
                                     compute_deformation_gradient_,
                                     "Compute and analyse the deformation gradient.")
                         .add_option("quadrature_order", quadrature_order, "Specify the quadrature order.")
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
            }

        private:
            FunctionSpace space_;
            Field_t field_;

            bool valid_{false};
            bool verbose_{false};
            bool compute_deformation_gradient_{false};
            int quadrature_order{0};
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_GRADIENT_APP_HPP

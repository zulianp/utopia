
#include "utopia_stk_intrepid2_Transport.hpp"
#include "utopia_intrepid2_Field.hpp"
#include "utopia_intrepid2_Gradient.hpp"
#include "utopia_intrepid2_Mass.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"
#include "utopia_intrepid2_Transport.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {
    namespace stk {

        class Transport::Impl {
        public:
            using Transport2 = utopia::Transport<2, Intrepid2FE::DynRankView>;
            using Transport3 = utopia::Transport<3, Intrepid2FE::DynRankView>;
            using DiffusionFunction = utopia::intrepid2::SubdomainValue<Scalar>;
            using Intrepid2Assembler = utopia::intrepid2::FEAssembler<Scalar>;

            std::shared_ptr<Field> field;
            std::shared_ptr<DiffusionFunction> diffusion_function;
            Scalar coeff{1.0};
            bool stabilize_transport{false};
            bool verbose{false};
            bool print_field{false};
        };

        Transport::Transport(const std::shared_ptr<FE> &fe) : Super(fe), impl_(utopia::make_unique<Impl>()) {}

        Transport::~Transport() = default;

        void Transport::read(Input &in) {
            Super::read(in);

            auto env = this->environment();

            if (env) {
                std::string pressure_field;
                in.require("pressure_field", pressure_field);

                auto p = env->find_field(*this->space(), pressure_field);

                if (!p) {
                    utopia::err() << "pressure_field not found in environment:";
                    env->describe(utopia::err().stream());
                }

                assert(p);

                if (p) {
                    impl_->field = p;
                } else {
                    std::string vector_field;
                    in.get("vector_field", vector_field);
                    auto v = env->find_field(*this->space(), vector_field);

                    if (v) {
                        impl_->field = v;
                    }
                }

            } else {
                std::string pressure_field;
                in.get("pressure_field", pressure_field);

                if (!pressure_field.empty()) {
                    assert(false);
                    Utopia::Abort("In order to retrive the pressure_field, The env must be defined!");
                }
            }

            in.get("coeff", impl_->coeff);
            in.get("stabilize_transport", impl_->stabilize_transport);
            in.get("verbose", impl_->verbose);
            in.get("print_field", impl_->print_field);

            const int spatial_dim = this->space()->mesh().spatial_dimension();

            intrepid2::Gradient<Scalar> g(this->fe_ptr());

            if (impl_->field->tensor_size() == 1) {
                assert(this->fe_ptr()->spatial_dimension() != 1);
                // If scalar field differentiate
                intrepid2::Field<Scalar> field(this->fe_ptr());
                convert_field(*impl_->field, field);
                g.init(field);
                g.scale(-impl_->coeff);

                in.get("diffusion_function", [this, &g](Input &node) {
                    impl_->diffusion_function = std::make_shared<utopia::intrepid2::SubdomainValue<Scalar>>(1.0);
                    // impl_->diffusion_function->read(node);
                    node.get("function", [this, &g](Input &inner_node) {
                        impl_->diffusion_function->read(inner_node);
                        g.scale(*impl_->diffusion_function);
                    });
                });

            } else {
                convert_field(*impl_->field, static_cast<intrepid2::Field<Scalar> &>(g));
            }

            if (impl_->print_field) {
                g.describe(utopia::out().stream());
            }

            switch (spatial_dim) {
                case 2: {
                    using Assemble2 = utopia::intrepid2::Assemble<Impl::Transport2>;
                    auto assembler = std::make_shared<Assemble2>(this->fe_ptr(), g.data());
                    assembler->read(in);
                    this->set_assembler(assembler);
                    break;
                }

                case 3: {
                    using Assemble3 = utopia::intrepid2::Assemble<Impl::Transport3>;
                    auto assembler = std::make_shared<Assemble3>(this->fe_ptr(), g.data());
                    assembler->read(in);
                    this->set_assembler(assembler);
                    break;
                }

                default: {
                    assert(false);
                    Utopia::Abort("utopia::stk::Transport: unsupported dimension " + std::to_string(spatial_dim) + "!");
                    break;
                }
            }

            if (impl_->verbose) {
                utopia::out() << "-----------------------------\n";
                utopia::out() << "Transport\n";
                if (impl_->field) {
                    utopia::out() << "Field:\t" << impl_->field->name() << '\n';
                }

                utopia::out() << "coeff:\t" << impl_->coeff << '\n';
                utopia::out() << "stabilize_transport:\t" << impl_->stabilize_transport << '\n';
                utopia::out() << "diffusion_function: ";
                // impl_->diffusion_function.describe(utopia::out().stream());
                utopia::out() << "print_field:\t" << impl_->print_field << '\n';
                utopia::out() << "-----------------------------\n";
            }
        }

    }  // namespace stk
}  // namespace utopia

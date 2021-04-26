
#include "utopia_stk_intrepid2_Transport.hpp"
#include "utopia_intrepid2_Field.hpp"
#include "utopia_intrepid2_Mass.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"
#include "utopia_intrepid2_Transport.hpp"
#include "utopia_stk_intrepid2.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {
    namespace stk {
        class StkIntrepid2Assembler::Impl {
        public:
            std::shared_ptr<Intrepid2FE> fe;
            std::shared_ptr<Intrepid2Assembler> assembler;

            inline bool empty() const { return static_cast<bool>(assembler); }
        };

        StkIntrepid2Assembler::~StkIntrepid2Assembler() = default;

        StkIntrepid2Assembler::StkIntrepid2Assembler() : impl_(utopia::make_unique<Impl>()) {}

        void StkIntrepid2Assembler::read(Input &) {}

        void StkIntrepid2Assembler::set_fe(const std::shared_ptr<Intrepid2FE> &fe) { impl_->fe = fe; }

        std::shared_ptr<StkIntrepid2Assembler::Intrepid2FE> StkIntrepid2Assembler::fe() { return impl_->fe; }

        void StkIntrepid2Assembler::ensure_fe(const int quadrature_order) {
            if (impl_->fe) {
                impl_->fe = std::make_shared<Intrepid2FE>();
                create_fe(*this->space(), *impl_->fe, quadrature_order);
            }
        }

        bool StkIntrepid2Assembler::assemble_element_tensors() {
            if (assembler()) {
                return assembler()->assemble();
            } else {
                return false;
            }
        }

        void StkIntrepid2Assembler::set_matrix_accumulator(const std::shared_ptr<TensorAccumulator> &accumulator) {
            assembler()->set_matrix_accumulator(accumulator);
        }

        std::shared_ptr<StkIntrepid2Assembler::TensorAccumulator> StkIntrepid2Assembler::matrix_accumulator() {
            return assembler()->matrix_accumulator();
        }

        void StkIntrepid2Assembler::set_assembler(const std::shared_ptr<Intrepid2Assembler> &assembler) {
            impl_->assembler = assembler;
        }

        std::shared_ptr<StkIntrepid2Assembler::Intrepid2Assembler> StkIntrepid2Assembler::assembler() {
            assert(impl_->assembler);
            return impl_->assembler;
        }

        bool StkIntrepid2Assembler::assemble(const Vector &x, Matrix &hessian, Vector &gradient) {
            if (!assemble_element_tensors()) {
                return false;
            }

            local_to_global(*space(), assembler()->matrix_accumulator()->data(), assembly_mode(), hessian);
            gradient = hessian * x;
            return true;
        }

        class Transport::Impl {
        public:
            using Transport2 = utopia::Transport<2, Intrepid2FE::DynRankView>;
            using Transport3 = utopia::Transport<3, Intrepid2FE::DynRankView>;
            using DiffusionFunction = utopia::intrepid2::SubdomainValue<Scalar>;

            std::shared_ptr<Field> field;
            std::shared_ptr<DiffusionFunction> diffusion_function;
            Scalar coeff{1.0};
            bool stabilize_transport{false};
            bool verbose{false};
            bool print_field{false};
        };

        Transport::Transport() : impl_(utopia::make_unique<Impl>()) {}

        Transport::~Transport() = default;

        void Transport::read(Input &in) {
            Super::read(in);
            this->ensure_fe(1);

            auto env = this->environment();

            if (env) {
                std::string field_field;
                in.get("pressure_field", field_field);
                auto p = env->find_field(*this->space(), field_field);

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
                std::string field_field;
                in.get("pressure_field", field_field);

                if (!field_field.empty()) {
                    assert(false);
                    Utopia::Abort("In order to retrive the field_field, The env must be defined!");
                }
            }

            in.get("coeff", impl_->coeff);
            in.get("stabilize_transport", impl_->stabilize_transport);
            in.get("verbose", impl_->verbose);
            in.get("print_field", impl_->print_field);

            const int spatial_dim = this->space()->mesh().spatial_dimension();

            intrepid2::Gradient<Scalar> g(this->fe());

            if (impl_->field->tensor_size() == 1) {
                assert(this->fe()->spatial_dimension() != 1);
                // If scalar field differentiate
                intrepid2::Field<Scalar> field(this->fe());
                convert_field(*impl_->field, field);
                g.init(field);
                g.scale(-impl_->coeff);

                in.get("diffusion_function", [this, &g](Input &in) {
                    impl_->diffusion_function = std::make_shared<utopia::intrepid2::SubdomainValue<Scalar>>(1.0);
                    impl_->diffusion_function->read(in);
                    g.scale(*impl_->diffusion_function);
                });

            } else {
                convert_field(*impl_->field, static_cast<intrepid2::Field<Scalar> &>(g));
            }

            switch (spatial_dim) {
                case 2: {
                    using Assemble2 = utopia::intrepid2::Assemble<Impl::Transport2>;
                    auto assembler = std::make_shared<Assemble2>(this->fe(), g.data());
                    assembler->read(in);
                    this->set_assembler(assembler);
                    break;
                }

                case 3: {
                    using Assemble3 = utopia::intrepid2::Assemble<Impl::Transport3>;
                    auto assembler = std::make_shared<Assemble3>(this->fe(), g.data());
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

        class Mass::Impl {
        public:
            using Mass = utopia::Mass<Scalar>;
            using Assemble = utopia::intrepid2::Assemble<Mass>;

            std::shared_ptr<Assemble> assembler;
            Impl() {}
        };

        Mass::~Mass() = default;
        Mass::Mass() : impl_(utopia::make_unique<Impl>()) {}

        void Mass::init() {}

        void Mass::ensure_assembler() {
            if (!impl_->assembler) {
                impl_->assembler = std::make_shared<Impl::Assemble>(this->fe());
            }
        }

        bool Mass::assemble_element_tensors() {
            if (!Super::assemble_element_tensors()) {
                return false;
            }

            // Check results
            return true;
        }

        void Mass::read(Input &in) {
            Super::read(in);
            impl_->assembler->read(in);
            init();
        }

    }  // namespace stk
}  // namespace utopia

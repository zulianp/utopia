
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
        class StkIntrepid2Assembler::Impl {
        public:
            std::shared_ptr<FunctionSpace> space;
            std::shared_ptr<Environment> environment;
        };

        void StkIntrepid2Assembler::set_environment(const std::shared_ptr<Environment> &env) {
            impl_->environment = env;
        }

        std::shared_ptr<StkIntrepid2Assembler::Environment> StkIntrepid2Assembler::environment() const {
            assert(impl_->environment);
            return impl_->environment;
        }

        void StkIntrepid2Assembler::set_space(const std::shared_ptr<FunctionSpace> &space) { impl_->space = space; }

        std::shared_ptr<FunctionSpace> StkIntrepid2Assembler::space() const { return impl_->space; }

        StkIntrepid2Assembler::~StkIntrepid2Assembler() = default;

        StkIntrepid2Assembler::StkIntrepid2Assembler(const std::shared_ptr<FE> &fe)
            : Super(fe), impl_(utopia::make_unique<Impl>()) {}

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

            std::shared_ptr<Intrepid2Assembler> assembler;
        };

        void Transport::set_matrix_accumulator(const std::shared_ptr<TensorAccumulator> &matrix_accumulator) {
            Super::set_matrix_accumulator(matrix_accumulator);
            assert(impl_->assembler);
            impl_->assembler->set_matrix_accumulator(matrix_accumulator);
        }

        void Transport::set_vector_accumulator(const std::shared_ptr<TensorAccumulator> &vector_accumulator) {
            Super::set_vector_accumulator(vector_accumulator);
            assert(impl_->assembler);
            impl_->assembler->set_vector_accumulator(vector_accumulator);
        }

        void Transport::set_scalar_accumulator(const std::shared_ptr<TensorAccumulator> &scalar_accumulator) {
            Super::set_scalar_accumulator(scalar_accumulator);
            assert(impl_->assembler);
            impl_->assembler->set_scalar_accumulator(scalar_accumulator);
        }

        void Transport::ensure_matrix_accumulator() {
            impl_->assembler->ensure_matrix_accumulator();
            Super::set_matrix_accumulator(impl_->assembler->matrix_accumulator());
        }

        void Transport::ensure_vector_accumulator() {
            impl_->assembler->ensure_vector_accumulator();
            Super::set_vector_accumulator(impl_->assembler->vector_accumulator());
        }

        void Transport::ensure_scalar_accumulator() {
            impl_->assembler->ensure_scalar_accumulator();
            Super::set_scalar_accumulator(impl_->assembler->scalar_accumulator());
        }

        Transport::Transport(const std::shared_ptr<FE> &fe) : Super(fe), impl_(utopia::make_unique<Impl>()) {}

        Transport::~Transport() = default;

        bool Transport::apply(const DynRankView &x, DynRankView &y) { return impl_->assembler->apply(x, y); }

        bool Transport::assemble_matrix() { return impl_->assembler->assemble_matrix(); }

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
                    impl_->assembler = assembler;
                    break;
                }

                case 3: {
                    using Assemble3 = utopia::intrepid2::Assemble<Impl::Transport3>;
                    auto assembler = std::make_shared<Assemble3>(this->fe_ptr(), g.data());
                    assembler->read(in);
                    impl_->assembler = assembler;
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


#include "utopia_stk_intrepid2_Transport.hpp"
#include "utopia_intrepid2_Field.hpp"
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

        void StkIntrepid2Assembler::read(Input &) {
            // if (!impl_->empty()) {
            //     assembler()->read(in);
            // } else {
            //     assert(false);
            // }
        }

        void StkIntrepid2Assembler::set_fe(const std::shared_ptr<Intrepid2FE> &fe) { impl_->fe = fe; }

        std::shared_ptr<StkIntrepid2Assembler::Intrepid2FE> StkIntrepid2Assembler::fe() { return impl_->fe; }

        void StkIntrepid2Assembler::ensure_fe(const int quadrature_order) {
            if (impl_->fe) {
                impl_->fe = std::make_shared<Intrepid2FE>();
                create_fe(*this->space(), *impl_->fe, quadrature_order);
            }
        }

        void StkIntrepid2Assembler::set_accumulator(const std::shared_ptr<TensorAccumulator> &accumulator) {
            assembler()->set_accumulator(accumulator);
        }

        std::shared_ptr<StkIntrepid2Assembler::TensorAccumulator> StkIntrepid2Assembler::accumulator() {
            return assembler()->accumulator();
        }

        void StkIntrepid2Assembler::set_assembler(const std::shared_ptr<Intrepid2Assembler> &assembler) {
            impl_->assembler = assembler;
        }

        std::shared_ptr<StkIntrepid2Assembler::Intrepid2Assembler> StkIntrepid2Assembler::assembler() {
            assert(impl_->assembler);
            return impl_->assembler;
        }

        class Transport::Impl {
        public:
            using Transport2 = utopia::Transport<2, Intrepid2FE::DynRankView>;
            using Transport3 = utopia::Transport<3, Intrepid2FE::DynRankView>;

            std::shared_ptr<Field> pressure;
            Scalar coeff{1.0};
            bool stabilize_transport{false};
            bool verbose{false};
            bool print_pressure{false};
        };

        Transport::Transport() : impl_(utopia::make_unique<Impl>()) {}

        Transport::~Transport() = default;

        bool Transport::assemble(const Vector &x, Matrix &hessian, Vector &gradient) { return false; }

        void Transport::read(Input &in) {
            Super::read(in);
            this->ensure_fe(1);

            auto env = this->environment();

            if (env) {
                std::string pressure_field;
                in.get("pressure_field", pressure_field);
                auto p = env->find_field(*this->space(), pressure_field);

                assert(p);

                if (p) {
                    impl_->pressure = p;
                }
            } else {
                std::string pressure_field;
                in.get("pressure_field", pressure_field);

                if (!pressure_field.empty()) {
                    assert(false);
                    Utopia::Abort("In order to retrive the pressure_field, The env must be defined!");
                }
            }

            // in.get("diffusion_function", impl_->diffusion_function);
            in.get("coeff", impl_->coeff);
            in.get("stabilize_transport", impl_->stabilize_transport);
            in.get("verbose", impl_->verbose);
            in.get("print_pressure", impl_->print_pressure);

            const int spatial_dim = this->space()->mesh().spatial_dimension();

            intrepid2::Field<Scalar> field(this->fe());
            convert_field(*impl_->pressure, field);

            intrepid2::Gradient<Scalar> g(this->fe());
            g.init(field);
            g.scale(-impl_->coeff);

            switch (spatial_dim) {
                case 2: {
                    using Assemble2 = utopia::intrepid2::Assemble<Impl::Transport2>;
                    auto assembler = std::make_shared<Assemble2>(g.data(), this->fe());
                    this->set_assembler(assembler);
                    break;
                }

                case 3: {
                    using Assemble3 = utopia::intrepid2::Assemble<Impl::Transport3>;
                    auto assembler = std::make_shared<Assemble3>(g.data(), this->fe());
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
                if (impl_->pressure) {
                    utopia::out() << "Pressure field:\t" << impl_->pressure->name() << '\n';
                }

                utopia::out() << "coeff:\t" << impl_->coeff << '\n';
                utopia::out() << "stabilize_transport:\t" << impl_->stabilize_transport << '\n';
                utopia::out() << "diffusion_function: ";
                // impl_->diffusion_function.describe(utopia::out().stream());
                utopia::out() << "print_pressure:\t" << impl_->print_pressure << '\n';
                utopia::out() << "-----------------------------\n";
            }
        }

    }  // namespace stk
}  // namespace utopia

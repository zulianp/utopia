#ifndef UTOPIA_KOKKOS_TRANSPORT_NEW_IMPL_HPP
#define UTOPIA_KOKKOS_TRANSPORT_NEW_IMPL_HPP

#include "utopia_kokkos_TransportNew.hpp"

#include "utopia_kokkos_Gradient.hpp"
#include "utopia_kokkos_Mass.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"
#include "utopia_kokkos_Transport.hpp"

#include "utopia_make_unique.hpp"

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE>
        class TransportNew<FunctionSpace, FE>::Impl {
        public:
            using DiffusionFunction = utopia::kokkos::SubdomainValue<FE>;
            using Field = utopia::kokkos::Field<FE>;
            using GradientField = utopia::kokkos::Gradient<FE>;
            using DynRankView = typename GradientField::DynRankView;

            using TransportOp2 = utopia::kokkos::kernels::
                TransportOp<2, Scalar, DynRankView, typename FE::Gradient, typename FE::Function, typename FE::Measure>;

            using TransportOp3 = utopia::kokkos::kernels::
                TransportOp<3, Scalar, DynRankView, typename FE::Gradient, typename FE::Function, typename FE::Measure>;

            std::shared_ptr<GradientField> vector_field;
            std::shared_ptr<DiffusionFunction> diffusion_function;
            Scalar coeff{1.0};
            bool stabilize_transport{false};
            bool verbose{false};
            bool print_field{false};
            int spatial_dimension{-1};
        };

        template <class FunctionSpace, class FE>
        TransportNew<FunctionSpace, FE>::TransportNew() : Super(), impl_(utopia::make_unique<Impl>()) {}

        template <class FunctionSpace, class FE>
        TransportNew<FunctionSpace, FE>::~TransportNew() = default;

        template <class FunctionSpace, class FE>
        void TransportNew<FunctionSpace, FE>::read(Input &in) {
            Super::read(in);
            // Model parameters
            in.get("coeff", impl_->coeff);
            in.get("stabilize_transport", impl_->stabilize_transport);
            in.get("verbose", impl_->verbose);
            in.get("print_field", impl_->print_field);

            // Vector field
            impl_->vector_field = std::make_shared<typename Impl::GradientField>(this->assembler()->fe_ptr());

            // Compute from pressure
            std::string pressure_field;
            in.require("pressure_field", pressure_field);
            auto p = this->field(pressure_field);

            if (p.empty()) {
                // Read directly from file
                std::string vector_field;
                in.get("vector_field", vector_field);
                auto v = this->field(vector_field);

                if (v.empty()) {
                    assert(false);
                    Utopia::Abort("Neither -- pressure_field -- nor -- vector_field -- are properly defined!");
                } else {
                    impl_->vector_field->data() = v[0]->data();
                }
            } else {
                impl_->vector_field->init(*p[0]);
                impl_->vector_field->scale(-impl_->coeff);

                in.get("diffusion_function", [this](Input &node) {
                    impl_->diffusion_function = std::make_shared<typename Impl::DiffusionFunction>(1.0);
                    node.get("function", [this](Input &inner_node) {
                        impl_->diffusion_function->read(inner_node);
                        impl_->vector_field->scale(*impl_->diffusion_function);
                    });
                });
            }

            auto discretization = this->assembler()->discretization();

            impl_->spatial_dimension = discretization->space()->mesh().spatial_dimension();

            if (impl_->print_field) {
                impl_->vector_field->describe(utopia::out().stream());
            }

            if (impl_->verbose) {
                utopia::out() << "-----------------------------\n";
                utopia::out() << "TransportNew\n";
                if (impl_->vector_field) {
                    utopia::out() << "Field:\t" << impl_->vector_field->name() << '\n';
                }

                utopia::out() << "coeff:\t" << impl_->coeff << '\n';
                utopia::out() << "stabilize_transport:\t" << impl_->stabilize_transport << '\n';
                utopia::out() << "diffusion_function: ";
                // impl_->diffusion_function.describe(utopia::out().stream());
                utopia::out() << "print_field:\t" << impl_->print_field << '\n';
                utopia::out() << "-----------------------------\n";
            }
        }

        template <class FunctionSpace, class FE>
        bool TransportNew<FunctionSpace, FE>::hessian_assemble(AssemblyMode mode) {
            UTOPIA_TRACE_REGION_BEGIN("TransportNew::hessian");

            auto &&assembler = this->assembler();
            assert(assembler);

            switch (impl_->spatial_dimension) {
                case 2: {
                    auto &&fe = assembler->fe();
                    typename Impl::TransportOp2 op2(impl_->vector_field->data(), fe.grad(), fe.fun(), fe.measure());
                    assembler->assemble_matrix_eij("TransportNew::hessian", mode, op2);
                    break;
                }
                case 3: {
                    auto &&fe = assembler->fe();
                    typename Impl::TransportOp3 op3(impl_->vector_field->data(), fe.grad(), fe.fun(), fe.measure());
                    assembler->assemble_matrix_eij("TransportNew::hessian", mode, op3);
                    break;
                }
                default: {
                    assert(false);
                    Utopia::Abort();
                }
            }

            UTOPIA_TRACE_REGION_END("TransportNew::hessian");
            return true;
        }

        template <class FunctionSpace, class FE>
        bool TransportNew<FunctionSpace, FE>::value_assemble(AssemblyMode mode) {
            assert(false);
            return false;
        }

        // Matrix free hessian application
        template <class FunctionSpace, class FE>
        bool TransportNew<FunctionSpace, FE>::apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) {
            UTOPIA_TRACE_REGION_BEGIN("TransportNew::apply");
            auto &&assembler = this->assembler();
            assert(assembler);

            switch (impl_->spatial_dimension) {
                case 2: {
                    auto &&fe = assembler->fe();
                    typename Impl::TransportOp2 op2(impl_->vector_field->data(), fe.grad(), fe.fun(), fe.measure());
                    assembler->assemble_apply_ei("TransportNew::apply", mode, op2, field);
                    break;
                }
                case 3: {
                    auto &&fe = assembler->fe();
                    typename Impl::TransportOp3 op3(impl_->vector_field->data(), fe.grad(), fe.fun(), fe.measure());
                    assembler->assemble_apply_ei("TransportNew::apply", mode, op3, field);
                    break;
                }
                default: {
                    assert(false);
                    Utopia::Abort();
                }
            }

            UTOPIA_TRACE_REGION_END("TransportNew::apply");
            return true;
        }

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_TRANSPORT_NEW_IMPL_HPP

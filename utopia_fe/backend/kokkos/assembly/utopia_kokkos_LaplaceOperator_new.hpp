#ifndef UTOPIA_KOKKOS_LAPLACE_OPERATOR_NEW_HPP
#define UTOPIA_KOKKOS_LAPLACE_OPERATOR_NEW_HPP

#include "utopia_fe_base.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_kokkos_LaplaceOp.hpp"
#include "utopia_kokkos_Material.hpp"

namespace utopia {

    namespace kokkos {

        template <class FunctionSpace, class FE_, class DiffusionCoefficient = typename FE_::Scalar>
        class LaplaceOperatorNew : public utopia::Material<FunctionSpace, FE_> {
        public:
            using FE = FE_;

            using Super = utopia::Material<FunctionSpace, FE_>;
            using Params = utopia::kokkos::LaplaceOp<FE, DiffusionCoefficient>;
            using Op = typename Params::UniformKernel;

            void read(Input &in) override { op_.read(in); }

            LaplaceOperatorNew(Params op = Params()) : Super(), op_(std::move(op)) {}

            inline int n_vars() const override { return 1; }
            inline std::string name() const override { return "LaplaceOperatorNew"; }

            inline bool has_hessian() const override { return true; }
            inline bool is_linear() const override { return true; }
            inline bool is_operator() const override { return true; };

            inline bool has_gradient() const override { return false; }
            inline bool has_value() const override { return false; }

            bool hessian_assemble(AssemblyMode mode) override {
                UTOPIA_TRACE_REGION_BEGIN("LaplaceOperatorNew::hessian");

                auto &&assembler = this->assembler();
                assert(assembler);

                assembler->assemble_matrix_eij(
                    "LaplaceOperatorNew::hessian", mode, op_.uniform_kernel(assembler->fe()));

                if (!op_.subdomain_value.empty) {
                    assert(false);
                }

                UTOPIA_TRACE_REGION_END("LaplaceOperatorNew::hessian");
                return true;
            }

            bool value_assemble(AssemblyMode mode) override {
                assert(false);
                return false;
            }

            // Matrix free hessian application
            bool apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) override {
                UTOPIA_TRACE_REGION_BEGIN("LaplaceOperatorNew::apply");
                auto &&assembler = this->assembler();
                assert(assembler);

                assembler->assemble_apply_ei(
                    "LaplaceOperatorNew::apply", mode, op_.uniform_kernel(assembler->fe()), field);

                UTOPIA_TRACE_REGION_END("LaplaceOperatorNew::apply");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_LAPLACE_OPERATOR_NEW_HPP

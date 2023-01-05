#ifndef UTOPIA_KOKKOS_FORCING_FUNCTION_NEW_HPP
#define UTOPIA_KOKKOS_FORCING_FUNCTION_NEW_HPP

#include "utopia_fe_base.hpp"

#include "utopia_Options.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_kokkos_MassOp.hpp"
#include "utopia_kokkos_Material.hpp"

namespace utopia {

    namespace kokkos {

        template <class FunctionSpace, class FE_, class DiffusionCoefficient = typename FE_::Scalar>
        class ForcingFunctionNew : public utopia::Material<FunctionSpace, FE_> {
        public:
            using Super = utopia::Material<FunctionSpace, FE_>;
            using FE = FE_;

            using Scalar = typename Traits<FE>::Scalar;

            using Measure = typename FE::Measure;
            using Fun = typename FE::Fun;

            class Params : public Configurable, public Describable {
            public:
                void read(Input &in) override {
                    Options()
                        .add_option("value", value, "Constant value for the force.")
                        .add_option("component", component, "Variable or component to which the force is applied.")
                        .add_option("verbose", verbose, "Verbose output.")
                        .add_option("n_components", n_components, "Number of components of the forcing function")
                        .add_option("density", density, "Domain density constant.")
                        .parse(in);
                }

                void describe(std::ostream &out = std::cout) const override {
                    out << "-----------------------------\n";
                    out << "ForcingFunction:\n";
                    out << "value:\t" << value << "\n";
                    out << "component:\t" << component << "\n";
                    out << "n_components:\t" << n_components << "\n";
                    out << "density:\t" << density << "\n";
                    out << "-----------------------------\n";
                }

                Params(const Scalar &value = Scalar(0.0)) : value(value) {}
                UTOPIA_FUNCTION Params(const Params &) = default;

                Scalar value;
                Scalar density{1.0};
                int n_components{1};
                int component{0};
                bool verbose{false};
            };

            void read(Input &in) override {
                Super::read(in);
                params_.read(in);
            }

            void describe(std::ostream &os = std::cout) const override {
                if (params_.verbose) {
                    if (!this->assembler()->space()->comm().rank()) {
                        params_.describe(os);
                    }
                }
            }

            ForcingFunctionNew(Params params = Params()) : Super(), params_(std::move(params)) {}

            inline int n_vars() const override { return 1; }
            inline std::string name() const override { return "ForcingFunction"; }

            inline bool has_hessian() const override { return false; }
            inline bool is_linear() const override { return true; }
            inline bool is_operator() const override { return false; };

            inline bool has_gradient() const override { return true; }
            inline bool has_value() const override { return false; }

            class FFOp : public TestOp {
            public:
                FFOp(const Scalar value, const Fun &fun, const Measure &measure, const int component)
                    : value(value), fun(fun), measure(measure) {
                    set_offset_test(component);
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i) const {
                    Scalar ret = 0.;

                    int n_qp = measure.extent(1);

                    for (int qp = 0; qp < n_qp; ++qp) {
                        auto dX = measure(cell, qp);

                        const Scalar f = fun(i, qp);

                        assert(f >= 0);
                        assert(f <= 1.0);

                        ret += -f * value * dX;
                    }

                    return ret;
                }

                Scalar value;
                Fun fun;
                Measure measure;
            };

            inline FFOp make_op() {
                auto &&assembler = this->assembler();
                auto fe = assembler->fe();
                return FFOp(params_.value * params_.density, fe.fun(), fe.measure(), params_.component);
            }

            bool gradient_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) override {
                UTOPIA_TRACE_REGION_BEGIN("ForcingFunctionNew::gradient");
                auto &&assembler = this->assembler();
                assert(assembler);

                assembler->assemble_vector_ei("ForcingFunctionNew::gradient", mode, make_op());

                UTOPIA_TRACE_REGION_END("ForcingFunctionNew::gradient");
                return true;
            }

            // NVCC_PRIVATE :
            Params params_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_FORCING_FUNCTION_NEW_HPP

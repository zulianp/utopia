
#ifndef UTOPIA_KOKKOS_FORCING_FUNCTION_HPP
#define UTOPIA_KOKKOS_FORCING_FUNCTION_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"

#include "utopia_Options.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, typename FE_, class Fun = typename FE_::Scalar>
        class ForcingFunction
            : public utopia::kokkos::FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using Scalar = typename Traits<Fun>::Scalar;
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>>;

            class Params : public Configurable {
            public:
                void read(Input &in) override {
                    Options()
                        .add_option("value", value, "Constant value for the force.")
                        .add_option("component", component, "Variable or component to which the force is applied.")
                        .add_option("verbose", verbose, "Verbose output.")
                        .add_option("n_components", n_components, "Number of components of the forcing function")
                        .add_option("density", density, "Domain density constant.")
                        .parse(in);

                    if (verbose) {
                        utopia::out() << "-----------------------------\n";
                        utopia::out() << "ForcingFunction:\n";
                        utopia::out() << "value:\t" << value << "\n";
                        utopia::out() << "component:\t" << component << "\n";
                        utopia::out() << "n_components:\t" << n_components << "\n";
                        utopia::out() << "density:\t" << density << "\n";
                        utopia::out() << "-----------------------------\n";
                    }
                }

                Params(const Fun &value = Fun(0.0)) : value(value) {}
                UTOPIA_FUNCTION Params(const Params &) = default;

                Fun value;
                Fun density{1.0};
                int n_components{1};
                int component{0};
                bool verbose{false};
            };

            ForcingFunction(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return op_.n_components; }

            inline bool is_matrix() const override { return false; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }

            inline std::string name() const override { return "ForcingFunction"; }

            bool assemble_vector() override {
                UTOPIA_TRACE_REGION_BEGIN("ForcingFunction::assemble");
                this->ensure_vector_accumulator();

                auto &fe = this->fe();

                // When integrating in subsets in parallel
                // we might get empty fe data
                if (!fe.empty()) {
                    auto ev = this->vector_data();

                    const int n_components = op_.n_components;
                    const int component = op_.component;
                    const int n_qp = fe.n_quad_points();

                    auto value = op_.value * op_.density;
                    auto fun = fe.fun();
                    auto measure = fe.measure();

                    this->loop_cell_test(
                        "ForcingFunction::assemble", KOKKOS_LAMBDA(const int &cell, const int &i) {
                            auto offset = i * n_components + component;

                            for (int qp = 0; qp < n_qp; ++qp) {
                                auto dX = measure(cell, qp);

                                const Scalar f = fun(i, qp);

                                assert(f >= 0);
                                assert(f <= 1.0);

                                assert(cell < int(ev.extent(0)));
                                assert(offset < int(ev.extent(1)));

                                ev(cell, offset) += -f * value * dX;
                            }
                        });

                    // if (op_.verbose) {
                    //     utopia::out() << "ForcingFunction: " << this->vector_accumulator()->sum() << '\n';
                    //     // this->describe(utopia::out().stream());
                    //     // utopia::out() << "Accumulator:\n";
                    //     // this->vector_accumulator()->describe(utopia::out().stream());
                    // }
                }

                UTOPIA_TRACE_REGION_END("ForcingFunction::assemble");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_FORCING_FUNCTION_HPP

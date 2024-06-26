
// #ifndef UTOPIA_KOKKOS_SYMBOLIC_FORCING_FUNCTION_HPP
// #define UTOPIA_KOKKOS_SYMBOLIC_FORCING_FUNCTION_HPP

// #include "utopia_kokkos_FE.hpp"
// #include "utopia_kokkos_FEAssembler.hpp"

// #include "utopia_Views.hpp"

// namespace utopia {
//     namespace kokkos {

//         template <typename FE_>
//         class SymbolicForcingFunction : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
//         public:
//             using Scalar = typename Traits<Fun>::Scalar;
//             using FE = FE_;
//             using SizeType = typename FE::SizeType;
//             using DynRankView = typename FE::DynRankView;
//             using FunctionSpaceTools = typename FE::FunctionSpaceTools;
//             using Fun = utopia::SymbolicFunction;

//             using ExecutionSpace = typename FE::ExecutionSpace;
//             using Super = utopia::kokkos::FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>>;

//             class Params : public Configurable {
//             public:
//                 void read(Input &in) override {
//                     in.get("value", fun);
//                     in.get("component", component);
//                     in.get("verbose", verbose);
//                     in.get("n_components", n_components);
//                 }

//                 Params(const Fun &fun = Fun("0.0")) : fun(fun) {}
//                 UTOPIA_FUNCTION Params(const Params &) = default;

//                 Fun fun;
//                 int n_components{1};
//                 int component{0};
//                 bool verbose{false};
//             };

//             SymbolicForcingFunction(const std::shared_ptr<FE> &fe, Params op = Params())
//                 : Super(fe), op_(std::move(op)) {}

//             inline int n_vars() const override { return op_.n_components; }

//             inline bool is_matrix() const override { return false; }
//             inline bool is_vector() const override { return true; }
//             inline bool is_scalar() const override { return false; }

//             inline std::string name() const override { return "SymbolicForcingFunction"; }

//             bool assemble_vector() override {
//                 UTOPIA_TRACE_REGION_BEGIN("Assemble<SymbolicForcingFunction>::assemble");
//                 this->ensure_vector_accumulator();

//                 auto &fe = this->fe();
//                 auto ev = this->vector_data();

//                 const int n_components = op_.n_components;
//                 const int component = op_.component;
//                 const int n_qp = fe.n_quad_points();

//                 auto value = op_.value;
//                 auto fun = fe.fun();
//                 auto measure = fe.measure();

//                 this->loop_cell_test(
//                     "Assemble<SymbolicForcingFunction>::assemble", KOKKOS_LAMBDA(const int &cell, const int &i) {
//                         auto offset = i * n_components + component;

//                         for (int qp = 0; qp < n_qp; ++qp) {
//                             auto dX = measure(cell, qp);

//                             const Scalar f = fun(i, qp);

//                             assert(f >= 0);
//                             assert(f <= 1.0);

//                             ev(cell, offset) += -f * value * dX;
//                         }
//                     });

//                 if (op_.verbose) {
//                     utopia::out() << "SymbolicForcingFunction: " << this->vector_accumulator()->sum() << '\n';
//                     // this->describe(utopia::out().stream());
//                     // utopia::out() << "Accumulator:\n";
//                     this->vector_accumulator()->describe(utopia::out().stream());
//                 }

//                 UTOPIA_TRACE_REGION_END("Assemble<SymbolicForcingFunction>::assemble");
//                 return true;
//             }

//             // NVCC_PRIVATE :
//             Params op_;
//         };
//     }  // namespace kokkos
// }  // namespace utopia

// #endif  // UTOPIA_KOKKOS_SYMBOLIC_FORCING_FUNCTION_HPP

// #ifndef UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP
// #define UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP

// #include "utopia_intrepid2_FE.hpp"
// #include "utopia_intrepid2_FEAssembler.hpp"
// #include "utopia_kokkos_LaplaceOp.hpp"

// #include "utopia_Views.hpp"

// namespace utopia {
//     template <int Dim, class Coefficient>
//     class VectorLaplaceOperator : public Configurable {
//     public:
//         void read(Input &in) override {
//             in.get("coeff", coeff);
//             in.get("verbose", verbose);
//         }

//         VectorLaplaceOperator(const Coefficient &coeff = Coefficient(1.0)) : coeff(coeff) {}
//         Coefficient coeff;
//         bool verbose{false};
//     };

//     namespace intrepid2 {

//         template <int Dim, typename Coefficient, typename Scalar>
//         class Assemble<VectorLaplaceOperator<Dim, Coefficient>, Scalar> : public FEAssembler<Scalar> {
//         public:
//             using FE = utopia::intrepid2::FE<Scalar>;
//             using SizeType = typename FE::SizeType;
//             using DynRankView = typename FE::DynRankView;
//             using FunctionSpaceTools = typename FE::FunctionSpaceTools;
//             using UserOp = utopia::VectorLaplaceOperator<Dim, Coefficient>;
//             using ExecutionSpace = typename FE::ExecutionSpace;
//             using Super = utopia::intrepid2::FEAssembler<Scalar>;

//             using Op = utopia::kokkos::kernels::
//                 VectorLaplaceOp<Dim, Scalar, Coefficient, typename FE::Gradient, typename FE::Measure>;

//             Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
//                 assert(Dim == fe->spatial_dimension());
//             }

//             inline int n_vars() const override { return Dim; }
//             inline std::string name() const override { return "VectorLaplaceOperator"; }

//             inline bool is_matrix() const override { return true; }
//             inline bool is_vector() const override { return true; }
//             inline bool is_scalar() const override { return false; }
//             bool is_operator() const override { return true; }

//             inline Op make_op() const { return Op(op_.coeff, this->fe().grad(), this->fe().measure()); }

//             bool apply(const DynRankView &x, DynRankView &y) override {
//                 UTOPIA_TRACE_REGION_BEGIN("Assemble<VectorLaplaceOperator>::apply");

//                 this->apply_vector_operator("Assemble<VectorLaplaceOperator>::apply", x, y, make_op());

//                 if (op_.verbose) {
//                     utopia::out() << "ForcingFunction: " << this->vector_accumulator()->sum() << '\n';
//                     // this->describe(utopia::out().stream());
//                     // utopia::out() << "Accumulator:\n";
//                     this->vector_accumulator()->describe(utopia::out().stream());
//                 }

//                 UTOPIA_TRACE_REGION_END("Assemble<VectorLaplaceOperator>::apply");
//                 return true;
//             }

//             bool assemble_matrix() override {
//                 UTOPIA_TRACE_REGION_BEGIN("Assemble<VectorLaplaceOperator>::assemble_matrix");

//                 this->ensure_matrix_accumulator();
//                 this->loop_cell_test_trial("Assemble<VectorLaplaceOperator>::assemble_matrix",
//                                            block_op_and_store_cell_ij(this->matrix_data(), make_op()));

//                 UTOPIA_TRACE_REGION_END("Assemble<VectorLaplaceOperator>::assemble_matrix");
//                 return true;
//             }

//             // NVCC_PRIVATE :
//             UserOp op_;
//         };
//     }  // namespace intrepid2
// }  // namespace utopia

// #endif  // UTOPIA_INTREPID2_VECTOR_LAPLACE_OPERATOR_HPP

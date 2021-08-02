// #ifndef UTOPIA_INTREPID2_TRANSPORT_HPP
// #define UTOPIA_INTREPID2_TRANSPORT_HPP

// #include "utopia_intrepid2_FE.hpp"
// #include "utopia_intrepid2_FEAssembler.hpp"
// #include "utopia_intrepid2_Gradient.hpp"

// #include "utopia_kokkos_TransportOp.hpp"

// #include "utopia_Views.hpp"

// namespace utopia {

//     template <int Dim, class Field>
//     class Transport : public Configurable {
//     public:
//         void read(Input &) override {}

//         Transport() = default;
//         Transport(const Field &vector_field) : vector_field(vector_field) {}
//         Field vector_field;
//     };

//     namespace intrepid2 {

//         template <int Dim, class Field>
//         class Assemble<Transport<Dim, Field>, typename Traits<Field>::Scalar>
//             : public utopia::intrepid2::FEAssembler<typename Traits<Field>::Scalar> {
//         public:
//             using Scalar = typename Traits<Field>::Scalar;
//             using FE = utopia::intrepid2::FE<Scalar>;
//             using SizeType = typename FE::SizeType;
//             using DynRankView = typename FE::DynRankView;
//             using UserOp = utopia::Transport<Dim, Field>;
//             using ExecutionSpace = typename FE::ExecutionSpace;
//             using Super = utopia::intrepid2::FEAssembler<Scalar>;

//             using Op = utopia::kokkos::kernels::
//                 TransportOp<Dim, Scalar, Field, typename FE::Gradient, typename FE::Function, typename FE::Measure>;

//             Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
//                 assert(Dim == fe->spatial_dimension());

// #ifndef NDEBUG
//                 if (op.vector_field.size() > 0) {
//                     assert(SizeType(op.vector_field.extent(0)) == fe->n_cells());
//                     assert(SizeType(op.vector_field.extent(1)) == fe->n_quad_points());
//                     assert(SizeType(op.vector_field.extent(2)) == fe->spatial_dimension());
//                 }
// #endif  // NDEBUG
//             }

//             inline int n_vars() const override { return 1; }

//             inline std::string name() const override { return "Transport"; }

//             inline bool is_matrix() const override { return true; }
//             inline bool is_vector() const override { return true; }
//             inline bool is_scalar() const override { return false; }
//             bool is_operator() const override { return true; }

//             inline Op make_op() {
//                 auto &fe = this->fe();
//                 return Op(op_.vector_field, fe.grad(), fe.fun(), fe.measure());
//             }

//             bool apply(const DynRankView &x, DynRankView &y) override {
//                 UTOPIA_TRACE_REGION_BEGIN("Assemble<Transport>::apply");

//                 this->apply_operator("Assemble<Transport>::apply", x, y, make_op());

//                 UTOPIA_TRACE_REGION_END("Assemble<Transport>::apply");
//                 return true;
//             }

//             bool assemble_matrix() override {
//                 UTOPIA_TRACE_REGION_BEGIN("Assemble<Transport>::assemble");

//                 this->ensure_matrix_accumulator();
//                 this->loop_cell_test_trial("Assemble<Transport>::assemble",
//                                            op_and_store_cell_ij(this->matrix_data(), make_op()));

//                 UTOPIA_TRACE_REGION_END("Assemble<Transport>::assemble");
//                 return true;
//             }

//             // NVCC_PRIVATE :
//             UserOp op_;
//         };
//     }  // namespace intrepid2

// }  // namespace utopia

// #endif  // UTOPIA_INTREPID2_TRANSPORT_HPP

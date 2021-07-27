#ifndef UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP
#define UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_FEAssembler.hpp"

#include "utopia_intrepid2_Strain.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    template <int Dim, class FirstLameParameter, class ShearModulus = FirstLameParameter>
    class LinearElasticity : public Configurable {
    public:
        void read(Input &in) override {
            StressStrainParameters<FirstLameParameter, ShearModulus> ssp;
            ssp.read(in);

            lambda = ssp.first_lame_parameter.get();
            mu = ssp.shear_modulus.get();
        }

        LinearElasticity(const FirstLameParameter &lambda = FirstLameParameter(1.0),
                         const ShearModulus &mu = FirstLameParameter(1.0))
            : lambda(lambda), mu(mu) {}

        FirstLameParameter lambda;
        ShearModulus mu;
    };

    namespace intrepid2 {

        template <int Dim, class FirstLameParameter, class ShearModulus, typename Scalar>
        class Assemble<LinearElasticity<Dim, FirstLameParameter, ShearModulus>, Scalar> : public FEAssembler<Scalar> {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using UserOp = utopia::LinearElasticity<Dim, FirstLameParameter, ShearModulus>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {
                assert(Dim == fe->spatial_dimension());
            }

            inline int n_vars() const override { return Dim; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            inline std::string name() const override { return "LinearElasticity"; }

            /// This is specialized for tensor-product finite elements
            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Op(const FirstLameParameter &lambda,
                                          const ShearModulus &mu,
                                          const DynRankView &grad,
                                          const DynRankView &measure)
                    : lambda(lambda), mux2(mu * 2), grad(grad), measure(measure), n_qp(measure.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar strain_inner(const int cell,
                                                           const int i,
                                                           const int j,
                                                           const int qp,
                                                           const int sub_i,
                                                           const int sub_j) const {
                    // StaticVector<Scalar, Dim> temp_i, temp_j;
                    // StaticMatrix<Scalar, Dim, Dim> strain_i;
                    // StaticMatrix<Scalar, Dim, Dim> strain_j;

                    // for (int d = 0; d < Dim; ++d) {
                    //     temp_i[d] = grad(cell, i, qp, d);
                    //     temp_j[d] = grad(cell, j, qp, d);
                    // }

                    // make_strain(sub_i, &temp_i[0], strain_i);
                    // make_strain(sub_j, &temp_j[0], strain_j);

                    // Scalar expected = inner(strain_i, strain_j);

                    auto ret = LinearizedStrain<Scalar, Dim>::inner(grad, cell, i, j, qp, sub_i, sub_j);
                    // assert(device::approxeq(ret, expected, 1e-10));
                    return ret;
                }

                UTOPIA_INLINE_FUNCTION Scalar strain_trace(const int cell,
                                                           const int i,
                                                           const int qp,
                                                           const int sub_i) const {
                    return LinearizedStrain<Scalar, Dim>::trace(grad, cell, i, qp, sub_i);
                }

                UTOPIA_INLINE_FUNCTION Scalar
                operator()(const int cell, const int i, const int j, const int sub_i, const int sub_j) const {
                    Scalar ret = 0.0;

                    for (int qp = 0; qp < n_qp; ++qp) {
                        const Scalar val = mux2 * strain_inner(cell, i, j, qp, sub_i, sub_j) +
                                           lambda * strain_trace(cell, i, qp, sub_i) * strain_trace(cell, j, qp, sub_j);

                        ret += val * measure(cell, qp);
                    }

                    return ret;
                }

                UTOPIA_INLINE_FUNCTION static constexpr int dim() { return Dim; }

                const FirstLameParameter lambda;
                const ShearModulus mux2;
                const DynRankView grad;
                const DynRankView measure;
                const int n_qp;
            };

            inline Op make_op() const { return Op(op_.lambda, op_.mu, this->fe().grad, this->fe().measure); }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LinearElasticity>::apply");

                this->apply_vector_operator("Assemble<LinearElasticity>::apply", x, y, make_op());

                UTOPIA_TRACE_REGION_END("Assemble<LinearElasticity>::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<LinearElasticity>::assemble_matrix");

                this->ensure_matrix_accumulator();
                this->loop_cell_test_trial("Assemble<LinearElasticity>::assemble_matrix",
                                           block_op_and_store_cell_ij(this->matrix_data(), make_op()));

                UTOPIA_TRACE_REGION_END("Assemble<LinearElasticity>::assemble_matrix");
                return true;
            }

            // NVCC_PRIVATE :
            UserOp op_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_LINEAR_ELASTICITY_HPP

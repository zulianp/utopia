#ifndef UTOPIA_INTREPID_2_MASS_HPP
#define UTOPIA_INTREPID_2_MASS_HPP

#include "utopia_intrepid2_LaplaceOperator.hpp"
#include "utopia_intrepid2_SubdomainFunction.hpp"

#include "utopia_Views.hpp"

namespace utopia {

    template <typename Fun>
    class Mass : public Configurable {
    public:
        using Scalar = typename Traits<Fun>::Scalar;
        // FIXME this type should be generalized to any backend
        using DensityFunction = utopia::intrepid2::SubdomainValue<Scalar>;

        void read(Input &in) override {
            in.get("density", density);
            in.get("n_components", n_components);
            in.get("lumped", lumped);
            in.get("verbose", verbose);
            in.get("expected_volume", expected_volume);
            in.get("expected_volume_tol", expected_volume_tol);

            in.get("density_function", [this](Input &in) {
                density_function = std::make_shared<DensityFunction>(1.0);
                density_function->read(in);
            });

            if (verbose) {
                utopia::out() << "-----------------------------\n";
                utopia::out() << "Mass\n";
                utopia::out() << "lumped:\t" << lumped << '\n';
                utopia::out() << "density:\t" << density << '\n';
                utopia::out() << "expected_volume:\t" << expected_volume << '\n';
                utopia::out() << "expected_volume_tol:\t" << expected_volume_tol << '\n';
                utopia::out() << "-----------------------------\n";
            }
        }

        Mass(const Fun &density = Fun(1.0)) : density(density) {}
        UTOPIA_FUNCTION Mass(const Mass &) = default;

        Fun density;
        int n_components{1};
        bool lumped{false};
        std::shared_ptr<DensityFunction> density_function;

        // Testing an printing
        bool verbose{false};
        Scalar expected_volume{-1};
        Scalar expected_volume_tol{1e-6};
    };

    namespace intrepid2 {

        template <typename Fun>
        class Assemble<Mass<Fun>, typename Traits<Fun>::Scalar> : public FEAssembler<typename Traits<Fun>::Scalar> {
        public:
            using Scalar = typename Traits<Fun>::Scalar;
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using UserOp = utopia::Mass<Fun>;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::intrepid2::FEAssembler<Scalar>;

            Assemble(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return op_.n_components; }
            int rank() const override { return 2; }
            inline std::string name() const override { return "Mass"; }

            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Op(const Fun &density,
                                          const DynRankView &fun,
                                          const DynRankView &measure,
                                          const int n_components)
                    : density(density),
                      fun(fun),
                      measure(measure),
                      n_components(n_components),
                      n_qp(measure.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i, const int &j) const {
                    auto offset_i = i * n_components;
                    auto offset_j = j * n_components;

                    Scalar ret = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        auto dX = measure(cell, qp);
                        ret += fun(offset_i, qp) * fun(offset_j, qp) * density * dX;
                    }

                    return ret;
                }

                const Fun density;
                const DynRankView fun;
                const DynRankView measure;
                const int n_components;
                const int n_qp;
            };

            template <class WrappedOp>
            class LumpedOp {
            public:
                UTOPIA_INLINE_FUNCTION LumpedOp(WrappedOp op) : op_(op), num_fields(op.fun.extent(1)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i) const {
                    Scalar ret = 0.0;
                    for (int j = 0; j < num_fields; ++j) {
                        ret += op_(cell, i, j);
                    }
                    return ret;
                }

                WrappedOp op_;
                const int num_fields;
            };

            class Lump {
            public:
                UTOPIA_INLINE_FUNCTION Lump(const DynRankView &element_matrices)
                    : element_matrices(element_matrices), num_fields(element_matrices.extent(2)) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int i, const int j) const {
                    if (i == j) {
                        return compute(cell, i);
                    } else {
                        return 0.0;
                    }
                }

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int i) const {
                    auto d = compute(cell, i);

                    for (int j = 0; j < num_fields; ++j) {
                        element_matrices(cell, i, j) = 0.0;
                    }

                    element_matrices(cell, i, i) = d;
                }

                UTOPIA_INLINE_FUNCTION Scalar compute(const int cell, const int i) const {
                    Scalar ret = 0.0;

                    for (int j = 0; j < num_fields; ++j) {
                        ret += element_matrices(cell, i, j);
                    }

                    return ret;
                }

                DynRankView element_matrices;
                const int num_fields;
            };

            inline Op make_op() {
                auto &fe = this->fe();
                return Op(op_.density, fe.fun, fe.measure, op_.n_components);
            }

            inline LumpedOp<Op> make_lumped_op() { return LumpedOp<Op>(make_op()); }

            bool apply(const DynRankView &x, DynRankView &y) override {
                if (op_.lumped) {
                    this->apply_diagonal_operator("Assemble<Mass>::apply::lumped", x, y, make_lumped_op());
                } else {
                    this->apply_operator("Assemble<Mass>::apply", x, y, make_op());
                }

                if (op_.density_function) {
                    this->scale_vector(*op_.density_function, y);
                }

                return true;
            }

            bool assemble() override {
                assert(op_.n_components == 1 && "IMPLEMENT ME");

                this->ensure_mat_accumulator();
                this->loop_cell_test_trial("Assemble<Mass>::assemble", op_and_store_cell_ij(this->data(), make_op()));

                if (op_.lumped) {
                    this->loop_cell_test("Assemble<Mass>::assemble::lump", Lump(this->data()));
                }

                if (op_.density_function) {
                    this->scale(*op_.density_function);
                }

                return true;
            }

            inline UserOp &user_op() { return op_; }
            inline const UserOp &user_op() const { return op_; }

            // NVCC_PRIVATE :
            UserOp op_;
        };
    }  // namespace intrepid2
}  // namespace utopia

#endif  // UTOPIA_INTREPID_2_MASS_HPP

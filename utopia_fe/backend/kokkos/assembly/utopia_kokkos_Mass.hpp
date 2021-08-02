#ifndef UTOPIA_KOKKOS_2_MASS_HPP
#define UTOPIA_KOKKOS_2_MASS_HPP

#include "utopia_kokkos_MassOp.hpp"

#include "utopia_kokkos_LaplaceOperator.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_Views.hpp"

namespace utopia {
    namespace kokkos {

        template <class FE_, typename Fun = typename FE_::Scalar>
        class Mass : public FEAssembler<FE_> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FE_>;

            using Op = utopia::kokkos::kernels::MassOp<Scalar, Fun, typename FE::Function, typename FE::Measure>;
            using LumpedOp = utopia::kokkos::kernels::LumpedOp<Op>;
            using Lump = utopia::kokkos::kernels::Lump<DynRankView>;

            class UserOp : public Configurable {
            public:
                using Scalar = typename Traits<Fun>::Scalar;
                // FIXME this type should be generalized to any backend
                // using DensityFunction = utopia::kokkos::SubdomainValue<Scalar>;
                using DensityFunction = utopia::kokkos::SubdomainValue<FE>;

                void read(Input &in) override {
                    in.get("density", density);
                    in.get("n_components", n_components);
                    in.get("lumped", lumped);
                    in.get("verbose", verbose);
                    in.get("expected_volume", expected_volume);
                    in.get("expected_volume_tol", expected_volume_tol);

                    in.get("density_function", [this](Input &node) {
                        density_function = std::make_shared<DensityFunction>(1.0);
                        node.get("function", [this](Input &inner_node) { density_function->read(inner_node); });
                    });

                    if (verbose) {
                        utopia::out() << "-----------------------------\n";
                        utopia::out() << "Mass\n";
                        utopia::out() << "lumped:\t" << lumped << '\n';
                        utopia::out() << "density:\t" << density << '\n';

                        if (expected_volume > 0) {
                            utopia::out() << "expected_volume:\t" << expected_volume << '\n';
                            utopia::out() << "expected_volume_tol:\t" << expected_volume_tol << '\n';
                        }

                        utopia::out() << "-----------------------------\n";
                    }
                }

                UserOp(const Fun &density = Fun(1.0)) : density(density) {}
                UTOPIA_FUNCTION UserOp(const UserOp &) = default;

                Fun density;
                int n_components{1};
                bool lumped{false};
                std::shared_ptr<DensityFunction> density_function;

                // Testing an printing
                bool verbose{false};
                Scalar expected_volume{-1};
                Scalar expected_volume_tol{1e-6};
            };

            Mass(const std::shared_ptr<FE> &fe, UserOp op = UserOp()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return op_.n_components; }
            inline std::string name() const override { return "Mass"; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            bool is_operator() const override { return true; }

            inline Op make_op() {
                auto &fe = this->fe();
                return Op(op_.density, fe.fun(), fe.measure(), op_.n_components);
            }

            inline LumpedOp make_lumped_op() { return LumpedOp(make_op()); }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<Mass>::apply");

                if (op_.lumped) {
                    if (op_.n_components == 1) {
                        this->apply_diagonal_operator("Assemble<Mass>::apply::lumped", x, y, make_lumped_op());
                    } else {
                        assert(false && "IMPLEMENT ME");
                    }
                } else {
                    if (op_.n_components == 1) {
                        this->apply_operator("Assemble<Mass>::apply", x, y, make_op());
                    } else {
                        this->apply_vector_operator("Assemble<Mass>::apply", x, y, make_op());
                    }
                }

                if (op_.density_function) {
                    this->scale_vector(*op_.density_function, y);
                }

                UTOPIA_TRACE_REGION_END("Assemble<Mass>::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<Mass>::assemble");

                this->ensure_matrix_accumulator();

                // assert(op_.n_components == 1 && "IMPLEMENT ME");

                if (op_.n_components == 1) {
                    this->loop_cell_test_trial("Assemble<Mass>::assemble_matrix",
                                               op_and_store_cell_ij(this->matrix_data(), make_op()));
                } else {
                    this->loop_cell_test_trial("Assemble<Mass>::assemble_matrix",
                                               block_op_and_store_cell_ij(this->matrix_data(), make_op()));
                }

                if (op_.lumped) {
                    this->loop_cell_test("Assemble<Mass>::assemble::lump", Lump(this->matrix_data()));
                }

                if (op_.density_function) {
                    auto data = this->matrix_data();
                    this->scale_matrix(*op_.density_function, data);
                }

                assert(check_volume());

                UTOPIA_TRACE_REGION_END("Assemble<Mass>::assemble");
                return true;
            }

            inline UserOp &user_op() { return op_; }
            inline const UserOp &user_op() const { return op_; }

            bool check_volume() const {
                if (op_.expected_volume <= 0) {
                    // Ignore
                    return true;
                }

                Scalar actual_volume = this->matrix_accumulator()->sum();

                if (!utopia::approxeq(actual_volume, op_.expected_volume, op_.expected_volume_tol)) {
                    utopia::err() << "Mass[Warning] expected volume " << op_.expected_volume << " have "
                                  << actual_volume << "instead!";
                    assert(false);
                    return false;
                }

                return true;
            }

            // NVCC_PRIVATE :
            UserOp op_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_2_MASS_HPP

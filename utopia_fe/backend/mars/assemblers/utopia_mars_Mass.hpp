#ifndef UTOPIA_MARS_MASS_HPP
#define UTOPIA_MARS_MASS_HPP

#include "mars_quad4.hpp"
#include "utopia_kokkos_MassOp.hpp"
#include "utopia_mars_ConcreteFEAssembler.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, typename... Args>
        class Mass : public ConcreteFEAssembler<DMesh, Args...> {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;
            using Scalar = Traits<mars::FunctionSpace>::Scalar;
            using Super = utopia::mars::ConcreteFEAssembler<DMesh, Args...>;
            static constexpr int Dim = DMesh::Dim;

            // Correct fe type is managed in ConcreteFEAssembler
            using FE = typename Super::FE;
            using Op = utopia::kokkos::kernels::MassOp<Scalar, Scalar, typename FE::Function, typename FE::Measure>;
            using LumpedOp = utopia::kokkos::kernels::LumpedOp<Op>;

            bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override {
                if (!assemble(hessian)) {
                    return false;
                }

                assemble(x, gradient);
                return true;
            }

            bool assemble(const Vector &, Matrix &hessian) override { return assemble(hessian); }

            bool assemble(const Vector &x, Vector &vec) override {
                this->ensure_fe();

                auto &fe = this->fe();

                // This is for linear materials since gradient := A x - b
                // non linear materials will have to implement a separate gradient Op

                Op op(density, fe.fun(), fe.measure(), n_components);

                if (lumped) {
                    Utopia::Abort("IMPLEMENT ME");
                    // return Super::block_op_apply_diagonal(LumpedOp(op), x, vec);
                    return false;
                } else {
                    return Super::block_op_apply(op, x, vec);
                }
            }

            bool assemble(Matrix &hessian) override {
                this->ensure_fe();

                auto &fe = this->fe();
                Op op(density, fe.fun(), fe.measure(), n_components);

                if (lumped) {
                    Utopia::Abort("IMPLEMENT ME");
                    // return Super::block_op_assemble_matrix_diagonal(LumpedOp(op), hessian);
                    return false;
                } else {
                    return Super::block_op_assemble_matrix(op, hessian);
                }
            }

            void read(Input &in) override {
                in.get("density", density);
                in.get("n_components", n_components);
                in.get("lumped", lumped);
                in.get("verbose", verbose);
                in.get("expected_volume", expected_volume);
                in.get("expected_volume_tol", expected_volume_tol);

                // FIXME
                // in.get("density_function", [this](Input &node) {
                //     density_function = std::make_shared<DensityFunction>(1.0);
                //     node.get("function", [this](Input &inner_node) { density_function->read(inner_node); });
                // });

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

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &env) override {
                Super::set_environment(env);
                // DO extra things ...
            }

            inline bool is_linear() const override { return true; }

            inline int n_vars() const { return n_components; }
            inline std::string name() const override { return "Mass"; }

        private:
            Scalar density;
            int n_components{1};
            bool lumped{false};

            // FIXME
            // std::shared_ptr<DensityFunction> density_function;

            // Testing an printing
            bool verbose{false};
            Scalar expected_volume{-1};
            Scalar expected_volume_tol{1e-6};
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_MASS_HPP

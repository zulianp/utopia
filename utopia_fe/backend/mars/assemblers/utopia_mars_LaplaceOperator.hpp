#ifndef UTOPIA_MARS_LAPLACE_OPERATOR_HPP
#define UTOPIA_MARS_LAPLACE_OPERATOR_HPP

#include "mars_quad4.hpp"
#include "utopia_mars_ConcreteFEAssembler.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, typename... Args>
        class LaplaceOperator : public ConcreteFEAssembler<DMesh, Args...> {
        public:
            using Matrix = Traits<mars::FunctionSpace>::Matrix;
            using Vector = Traits<mars::FunctionSpace>::Vector;
            using Scalar = Traits<mars::FunctionSpace>::Scalar;
            using Super = utopia::mars::ConcreteFEAssembler<DMesh, Args...>;

            using DofHandler = typename Super::DofHandler;
            using FEDofMap = typename Super::FEDofMap;
            using SPattern = typename Super::SPattern;

            static const ::mars::Integer Type = SPattern::DofHandler::ElemType;

            using UniformFE = typename utopia::mars::FETypeSelect<Scalar, DMesh::Dim>::Type;

            bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override {
                if (!assemble(hessian)) {
                    return false;
                }

                gradient = hessian * x;
                return true;
            }

            bool assemble(const Vector &, Matrix &hessian) override { return assemble(hessian); }

            bool assemble(const Vector &x, Vector &vec) override {
                assert(false && "IMPLEMENT ME");
                return false;
            }

            bool assemble(Matrix &hessian) override {
                init();

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                auto measure = fe_->measure;
                auto grad = fe_->grad;

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < UniformFE::NFun; i++) {
                        const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, i);
                        if (!dof_handler.is_owned(local_dof_i)) continue;

                        for (int j = 0; j < UniformFE::NFun; j++) {
                            const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, j);
                            // for each dof get the local number

                            Scalar val = 0.0;
                            for (int k = 0; k < UniformFE::NQPoints; ++k) {
                                Scalar dot_grads = 0.0;

                                for (int d = 0; d < DMesh::Dim; ++d) {
                                    assert(grad(i, k, d) == grad(i, k, d));
                                    assert(grad(j, k, d) == grad(j, k, d));

                                    auto g_i = grad(i, k, d);
                                    auto g_j = grad(j, k, d);

                                    dot_grads += g_i * g_j;
                                    assert(dot_grads == dot_grads);
                                }

                                Scalar dX = measure(k);
                                assert(dX == dX);

                                val += dot_grads * dX;
                            }

                            assert(val == val);
                            sp.atomic_add_value(local_dof_i, local_dof_j, val);
                        }
                    }
                });

                // // FIXME Bad it should not do this
                auto mat_impl =
                    Teuchos::rcp(new Matrix::CrsMatrixType(this->handler()->get_sparsity_pattern().get_matrix(),
                                                           hessian.raw_type()->getRowMap(),
                                                           hessian.raw_type()->getColMap()));
                hessian.wrap(mat_impl, false);
                return true;
            }

            void read(Input &in) override { in.get("coeff", coeff_); }

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &) override {
                // assert(false);
            }

            inline bool is_linear() const override { return true; }

        private:
            Scalar coeff_{1.0};
            std::unique_ptr<UniformFE> fe_;

            void init() {
                if (!fe_) {
                    using Coordinates = typename UniformFE::Coordinates;

                    Coordinates coords("coords");

                    auto handler = this->handler();
                    auto sp = handler->get_sparsity_pattern();
                    auto fe_dof_map = handler->get_fe_dof_map();
                    auto dof_handler = sp.get_dof_handler();

                    fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                        if (elem_index == 0) {
                            Scalar p[DMesh::Dim];

                            for (int i = 0; i < coords.extent(0); ++i) {
                                auto local_dof = fe_dof_map.get_elem_local_dof(elem_index, i);
                                dof_handler.template get_dof_coordinates_from_local<Type>(local_dof, p);

                                for (int d = 0; d < coords.extent(1); ++d) {
                                    coords(i, d) = p[d];
                                }
                            }
                        }
                    });

                    fe_ = utopia::make_unique<UniformFE>();
                    fe_->init(coords);
                    //     auto handler = this->handler();
                    //     fe_ = utopia::make_unique<FE>(handler->get_sparsity_pattern(), handler->get_fe_dof_map());
                    //     fe_->init();
                }
            }
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_LAPLACE_OPERATOR_HPP

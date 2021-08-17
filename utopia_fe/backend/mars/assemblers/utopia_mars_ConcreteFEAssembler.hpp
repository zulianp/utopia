#ifndef UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP
#define UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP

#include "utopia_mars_ForwardDeclarations.hpp"

#include "utopia_mars_FEAssembler.hpp"
#include "utopia_mars_FEHandler.hpp"

#include "utopia_mars_FE.hpp"
#include "utopia_mars_UniformFE.hpp"

#include "utopia_mars_FEBuilder.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, typename...>
        class ConcreteFEAssembler : public FEAssembler {
        public:
            using Scalar = typename Traits<FunctionSpace>::Scalar;

            using FEHandler = utopia::mars::FEHandler<DMesh, 1>;
            using DofHandler = typename FEHandler::DofHandler;
            using FEDofMap = typename FEHandler::FEDofMap;
            using SPattern = typename FEHandler::SPattern;

            static const ::mars::Integer Type = SPattern::DofHandler::ElemType;
            using UniformFE = utopia::kokkos::UniformFE<Scalar>;

            // using FE = utopia::mars::FE<Scalar, SPattern, FEDofMap>;
            using FE = UniformFE;

            virtual ~ConcreteFEAssembler() = default;

            inline void set_space(const std::shared_ptr<FunctionSpace> &space) override { space_ = space; }
            inline std::shared_ptr<FunctionSpace> space() const { return space_; }
            // inline std::shared_ptr<MarsMeshType> mars_mesh_ptr() { return space_->mesh().raw_type<MarsMeshType>(); }
            inline std::shared_ptr<FEHandler> handler() { return space_->raw_type<FEHandler>(); }

            template <class Op>
            bool scalar_op_assemble_matrix(Op op, Matrix &hessian) {
                ensure_fe();

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, i);
                        if (!dof_handler.is_owned(local_dof_i)) continue;

                        for (int j = 0; j < n_fun; j++) {
                            const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, j);
                            Scalar val = op(elem_index, i, j);

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

            template <class Op>
            bool scalar_op_apply(Op op, const Vector &x, Vector &y) {
                ensure_fe();

                if (empty(y)) {
                    y.zeros(layout(x));
                } else {
                    y.set(0.0);
                }

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                // FIXME Missing gobal to local for x!!!
                auto x_view = local_view_device(x).raw_type();
                auto y_view = local_view_device(y).raw_type();

                // kokkos_view_owned = local_view_device(x).raw_type(); // Rank 2 tensor N x 1
                // View x_local; // Rank 1 tensor
                // fe_dof_map.copy_owned_values(kokkos_view_owned, x_local);
                // fe_dof_map.gather_values(kokkos_view_owned, x_local);

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, i);
                        if (!dof_handler.is_owned(local_dof_i)) continue;

                        Scalar val = 0;

                        for (int j = 0; j < n_fun; j++) {
                            // FIXME What is the correct index map?
                            const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, j);
                            auto x_j = x_view(local_dof_j, 0);

                            val += op(elem_index, i, j) * x_j;

                            // val += op(elem_index, i, j, block_i, block_j) * x_local_j;
                        }

                        Kokkos::atomic_fetch_add(&y_view(local_dof_i, 0), val);
                    }
                });

                return true;
            }

            void read(Input &in) override {}

            void set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &) override {
                // assert(false);
            }

            inline bool is_linear() const override { return true; }

            inline UniformFE &fe() {
                assert(fe_);
                return *fe_;
            }

            void ensure_fe() { init(); }

        private:
            std::shared_ptr<FunctionSpace> space_;
            std::unique_ptr<UniformFE> fe_;

            void init() {
                if (!fe_) {
                    fe_ = utopia::make_unique<UniformFE>();
                    FEBuilder<FEHandler, utopia::kokkos::UniformFE<Scalar>> builder;
                    builder.build(*this->handler(), *fe_);
                }
            }
        };

    }  // namespace mars
}  // namespace utopia

#endif  // UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP

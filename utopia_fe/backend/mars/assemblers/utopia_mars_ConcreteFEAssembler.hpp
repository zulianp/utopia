#ifndef UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP
#define UTOPIA_MARS_CONCRETE_FE_ASSEMBLER_HPP

#include "utopia_mars_ForwardDeclarations.hpp"

#include "utopia_mars_FEAssembler.hpp"
#include "utopia_mars_FEHandler.hpp"

#include "utopia_mars_FE.hpp"
#include "utopia_mars_UniformFE.hpp"

#include "utopia_mars_FEBuilder.hpp"

#include "mars_distributed_base_data_management.hpp"

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
            using FEAssembler::assemble;

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

                if (!add_offsetted_scalar_op_to_matrix(0, 0, op, sp)) {
                    return false;
                }

                // FIXME Bad it should not do this
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

                auto y_view = local_view_device(y).raw_type();                                  // Rank 2 tensor N x 1
                ::mars::ViewVectorType<Scalar> x_local("x_local", dof_handler.get_dof_size());  // Rank 1 tensor
                collect_ghost_layer(x, x_local);

                return add_offsetted_scalar_op_to_vector(0, 0, op, x_local, y_view);
            }

            template <class Op, class TpetraVectorView>
            bool add_offsetted_scalar_op_to_vector(const int s_offset_i,
                                                   const int s_offset_j,
                                                   Op op,
                                                   const ::mars::ViewVectorType<Scalar> &x,
                                                   TpetraVectorView &y) {
                ensure_fe();

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        auto offset_i = dof_handler.compute_block_index(i, s_offset_i);
                        const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                        if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                        Scalar val = 0;
                        for (int j = 0; j < n_fun; j++) {
                            auto offset_j = dof_handler.compute_block_index(i, s_offset_j);
                            const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);
                            if (local_dof_j > -1) {
                                assert(local_dof_j < x.extent(0));
                                auto x_j = x(local_dof_j);
                                val += op(elem_index, i, j) * x_j;
                            }
                        }

                        const auto owned_dof_i = dof_handler.local_to_owned_index(local_dof_i);
                        Kokkos::atomic_fetch_add(&y(owned_dof_i, 0), val);
                    }
                });

                return true;
            }

            void collect_ghost_layer(const Vector &in, ::mars::ViewVectorType<Scalar> &out) {
                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                auto x_view = local_view_device(in).raw_type();  // Rank 2 tensor N x 1
                if (out.extent(0) != dof_handler.get_dof_size()) {
                    out = ::mars::ViewVectorType<Scalar>("x_local", dof_handler.get_dof_size());  // Rank 1 tensor
                }

                ::mars::ViewVectorType<Scalar> x_view_rank1(::Kokkos::subview(x_view, ::Kokkos::ALL, 0));
                ::mars::set_locally_owned_data(dof_handler, out, x_view_rank1);
                ::mars::gather_ghost_data(dof_handler, out);
            }

            ////////////////////////

            template <class Op>
            bool block_op_assemble_matrix(Op op, Matrix &hessian) {
                ensure_fe();

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                int block_size = op.dim();

                if (block_size != dof_handler.get_block()) {
                    utopia::err() << block_size << "==" << dof_handler.get_block() << "\n";
                }

                assert(block_size == dof_handler.get_block());

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        for (int sub_i = 0; sub_i < block_size; ++sub_i) {
                            auto offset_i = dof_handler.compute_block_index(i, sub_i);
                            const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                            if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                            for (int j = 0; j < n_fun; j++) {
                                for (int sub_j = 0; sub_j < block_size; ++sub_j) {
                                    auto offset_j = dof_handler.compute_block_index(j, sub_j);
                                    const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);
                                    if (local_dof_j > -1) {
                                        const Scalar val = op(elem_index, i, j, sub_i, sub_j);

                                        assert(val == val);
                                        sp.atomic_add_value(local_dof_i, local_dof_j, val);
                                    }
                                }
                            }
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

            /////////////////

            template <class Op>
            bool block_op_apply(Op op, const Vector &x, Vector &y) {
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

                auto x_view = local_view_device(x).raw_type();
                auto y_view = local_view_device(y).raw_type();

                // kokkos_view_owned = local_view_device(x).raw_type(); // Rank 2 tensor N x 1

                ::mars::ViewVectorType<Scalar> x_local("x_local", dof_handler.get_dof_size());  // Rank 1 tensor
                collect_ghost_layer(x, x_local);

                return add_offsetted_block_op_to_vector(0, 0, op, x_local, y_view);
            }

            template <class Op, class TpetraVectorView>
            bool add_offsetted_block_op_to_vector(const int b_offset_i,
                                                  const int b_offset_j,
                                                  Op op,
                                                  const ::mars::ViewVectorType<Scalar> &x,
                                                  TpetraVectorView &y) {
                ensure_fe();

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                int block_size = op.dim();
                assert(block_size <= dof_handler.get_block());

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        for (int sub_i = 0; sub_i < block_size; ++sub_i) {
                            auto offset_i = dof_handler.compute_block_index(i, sub_i + b_offset_i);
                            const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                            if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                            Scalar val = 0;
                            for (int j = 0; j < n_fun; j++) {
                                for (int sub_j = 0; sub_j < block_size; ++sub_j) {
                                    auto offset_j = dof_handler.compute_block_index(j, sub_j + b_offset_j);
                                    const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);
                                    if (local_dof_i > -1) {
                                        auto x_j = x(local_dof_j);
                                        val += op(elem_index, i, j, sub_i, sub_j) * x_j;
                                    }
                                }
                            }

                            const auto owned_dof_i = dof_handler.local_to_owned_index(local_dof_i);
                            Kokkos::atomic_fetch_add(&y(owned_dof_i, 0), val);
                        }
                    }
                });

                return true;
            }

            template <class Op, class TpetraVectorView>
            bool add_offsetted_block_x_scalar_op_to_vector(const int b_offset_i,
                                                           const int s_offset_j,
                                                           Op op,
                                                           const ::mars::ViewVectorType<Scalar> &x,
                                                           TpetraVectorView &y) {
                ensure_fe();

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                int block_size = op.dim();
                assert(block_size <= dof_handler.get_block());

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        for (int sub_i = 0; sub_i < block_size; ++sub_i) {
                            auto offset_i = dof_handler.compute_block_index(i, sub_i + b_offset_i);
                            const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                            if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                            Scalar val = 0;
                            for (int j = 0; j < n_fun; j++) {
                                // for (int sub_j = 0; sub_j < block_size; ++sub_j) {
                                auto offset_j = dof_handler.compute_block_index(j, s_offset_j);
                                const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);
                                if (local_dof_i > -1) {
                                    auto x_j = x(local_dof_j);
                                    val += op(elem_index, i, j, sub_i) * x_j;
                                }
                                // }
                            }

                            const auto owned_dof_i = dof_handler.local_to_owned_index(local_dof_i);
                            Kokkos::atomic_fetch_add(&y(owned_dof_i, 0), val);
                        }
                    }
                });

                return true;
            }

            template <class VecOp, class ScalarOp, class CoupleVecScalarOp>
            bool coupled_vector_scalar_op_assemble(const int vector_var,
                                                   const int scalar_var,
                                                   VecOp vec_op,
                                                   ScalarOp scalar_op,
                                                   CoupleVecScalarOp couple_vec_scalar_op,
                                                   Matrix &hessian) {
                ensure_fe();

                auto handler = this->handler();
                auto sp = handler->get_sparsity_pattern();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                int block_size = vec_op.dim();
                assert(block_size <= dof_handler.get_block());

                bool ok = add_offsetted_block_op_to_matrix(vector_var, vector_var, vec_op, sp) &&
                          add_offsetted_scalar_op_to_matrix(scalar_var, scalar_var, scalar_op, sp) &&
                          add_offsetted_block_x_scalar_op_to_matrix(vector_var, scalar_var, couple_vec_scalar_op, sp);

                if (!ok) return false;

                // // FIXME Bad it should not do this
                auto mat_impl =
                    Teuchos::rcp(new Matrix::CrsMatrixType(this->handler()->get_sparsity_pattern().get_matrix(),
                                                           hessian.raw_type()->getRowMap(),
                                                           hessian.raw_type()->getColMap()));
                hessian.wrap(mat_impl, false);
                return true;
            }

            template <class VecOp, class ScalarOp, class CoupleVecScalarOp>
            bool coupled_vector_scalar_op_apply(const int vector_var,
                                                const int scalar_var,
                                                VecOp vec_op,
                                                ScalarOp scalar_op,
                                                CoupleVecScalarOp couple_vec_scalar_op,
                                                const Vector &x,
                                                Vector &y) {
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

                // Why does it crash here???
                ::mars::ViewVectorType<Scalar> weird_x_local("weird_x_local",
                                                             dof_handler.get_dof_size());  // Rank 1 tensor
                collect_ghost_layer(x, weird_x_local);

                auto y_view = local_view_device(y).raw_type();

                int block_size = vec_op.dim();
                assert(block_size <= dof_handler.get_block());

                return add_offsetted_block_op_to_vector(vector_var, vector_var, vec_op, weird_x_local, y_view) &&
                       add_offsetted_block_x_scalar_op_to_vector(
                           vector_var, scalar_var, couple_vec_scalar_op, weird_x_local, y_view) &&
                       add_offsetted_scalar_op_to_vector(scalar_var, scalar_var, scalar_op, weird_x_local, y_view);
            }

            ////////////////////

            template <class Op>
            bool add_offsetted_scalar_op_to_matrix(const int s_offset_i, const int s_offset_j, Op op, SPattern &sp) {
                ensure_fe();

                auto handler = this->handler();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        auto offset_i = dof_handler.compute_block_index(i, s_offset_i);
                        const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                        if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                        for (int j = 0; j < n_fun; j++) {
                            auto offset_j = dof_handler.compute_block_index(i, s_offset_j);
                            const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);
                            if (local_dof_j > -1) {
                                Scalar val = op(elem_index, i, j);

                                assert(val == val);
                                sp.atomic_add_value(local_dof_i, local_dof_j, val);
                            }
                        }
                    }
                });

                return true;
            }

            ////////////////////

            template <class Op>
            bool add_offsetted_block_op_to_matrix(const int b_offset_i, const int b_offset_j, Op op, SPattern &sp) {
                ensure_fe();

                auto handler = this->handler();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                int block_size = op.dim();

                if (block_size > dof_handler.get_block()) {
                    utopia::err() << block_size << "<=" << dof_handler.get_block() << "\n";
                }

                assert(block_size <= dof_handler.get_block());

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        for (int sub_i = 0; sub_i < block_size; ++sub_i) {
                            auto offset_i = dof_handler.compute_block_index(i, b_offset_i + sub_i);
                            const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                            if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                            for (int j = 0; j < n_fun; j++) {
                                for (int sub_j = 0; sub_j < block_size; ++sub_j) {
                                    auto offset_j = dof_handler.compute_block_index(j, b_offset_j + sub_j);
                                    const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);
                                    if (local_dof_j > -1) {
                                        const Scalar val = op(elem_index, i, j, sub_i, sub_j);

                                        assert(val == val);
                                        sp.atomic_add_value(local_dof_i, local_dof_j, val);
                                    }
                                }
                            }
                        }
                    }
                });

                return true;
            }

            ////////////////////

            template <class Op>
            bool add_offsetted_block_x_scalar_op_to_matrix(const int b_offset_i,
                                                           const int s_offset_j,
                                                           Op op,
                                                           SPattern &sp) {
                ensure_fe();

                auto handler = this->handler();
                auto fe_dof_map = handler->get_fe_dof_map();
                auto dof_handler = sp.get_dof_handler();

                const int n_fun = fe_->n_shape_functions();
                const int n_qp = fe_->n_quad_points();

                int block_size = op.dim();

                if (block_size > dof_handler.get_block()) {
                    utopia::err() << block_size << "<=" << dof_handler.get_block() << "\n";
                }

                assert(block_size <= dof_handler.get_block());

                fe_dof_map.iterate(MARS_LAMBDA(const ::mars::Integer elem_index) {
                    for (int i = 0; i < n_fun; i++) {
                        for (int sub_i = 0; sub_i < block_size; ++sub_i) {
                            auto offset_i = dof_handler.compute_block_index(i, b_offset_i + sub_i);
                            const auto local_dof_i = fe_dof_map.get_elem_local_dof(elem_index, offset_i);
                            if (local_dof_i < 0 || !dof_handler.is_owned(local_dof_i)) continue;

                            for (int j = 0; j < n_fun; j++) {
                                auto offset_j = dof_handler.compute_block_index(j, s_offset_j);
                                const auto local_dof_j = fe_dof_map.get_elem_local_dof(elem_index, offset_j);
                                if (local_dof_j > -1) {
                                    const Scalar val = op(elem_index, i, j, sub_i);

                                    assert(val == val);
                                    sp.atomic_add_value(local_dof_i, local_dof_j, val);
                                }
                            }
                        }
                    }
                });

                return true;
            }

            void read(Input &in) override {}

            void set_environment(const std::shared_ptr<Environment> &) override {
                // assert(false);
            }

            inline bool is_linear() const override { return true; }

            inline UniformFE &fe() {
                assert(fe_);
                return *fe_;
            }

            void ensure_fe() { init(); }

            std::string name() const override {
                Utopia::Abort("IMPLEMENT ME");
                return "Undefined!";
            }

            bool apply(const Vector &, Vector &) override {
                Utopia::Abort("IMPLEMENT ME");
                return false;
            }

            bool assemble(Vector &) override {
                Utopia::Abort("IMPLEMENT ME");
                return false;
            }

            std::shared_ptr<Environment> environment() const override {
                Utopia::Abort("IMPLEMENT ME");
                return nullptr;
            }

            void set_time(const std::shared_ptr<SimulationTime> &time) override { Utopia::Abort("IMPLEMENT ME"); }

            bool assemble(const Vector &x, Vector &vec) override {
                assert(is_linear());

                if (empty(vec)) {
                    vec.zeros(layout(x));
                }

                bool ok = apply(x, vec);
                assert(ok);
                return ok;
            }

            bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override {
                if (!this->assemble(hessian)) {
                    return false;
                }

                return this->assemble(x, gradient);
            }

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

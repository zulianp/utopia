#ifndef UTOPIA_MARS_FEHANDLER_HPP
#define UTOPIA_MARS_FEHANDLER_HPP

#include "mars.hpp"
#include "mars_base.hpp"
#include "mars_context.hpp"
#include "mars_distributed_dof_management.hpp"
#include "mars_globals.hpp"
#include "utopia_mars_FunctionSpace.hpp"

#include "utopia_mars_Factory_impl.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, int Degree_ = 1>
        class FEHandler : public IFEHandler {
        public:
            static constexpr ::mars::Integer Degree = Degree_;
            static constexpr int Dim = DMesh::Dim;
            using DofHandler = ::mars::DofHandler<DMesh, Degree, 0>;
            // using DofHandler = ::mars::DofHandler<DMesh, Degree, 1>;
            using FEDofMap = ::mars::FEDofMap<DofHandler>;
            using SPattern =
                ::mars::SparsityPattern<Scalar, LocalSizeType, SizeType, DofHandler, MarsCrsMatrix::size_type>;

            static_assert(std::is_same<SizeType, Matrix::CrsMatrixType::global_ordinal_type>::value, "Weird!");

            auto new_crs_matrix() -> MarsCrsMatrix override { return sparsity_pattern->new_crs_matrix(); }

            void describe() const override {
                dof_handler->print_dofs();
                sparsity_pattern->print_sparsity_pattern();

                // ::mars::print_fe_dof_map(*dof_handler_impl, *fe_dof_map_impl);
            }

            void matrix_apply_constraints(Matrix &m, const Scalar diag_value, const std::string side) override {
                // BC set constrained rows to zero, except diagonal where you set diag_value
                auto sp = *sparsity_pattern;
                auto mat = m.raw_type()->getLocalMatrix();
                dof_handler->boundary_dof_iterate(
                    MARS_LAMBDA(const ::mars::Integer local_dof) {
                        sp.matrix_apply_constraints(local_dof, mat, diag_value);
                    },
                    side);
            }

            void vector_apply_constraints(Vector &v, const Scalar value, const std::string side) override {
                auto sp = *sparsity_pattern;
                auto vec = v.raw_type()->getLocalView<::mars::KokkosSpace>();
                // BC set values to constraint value (i.e., boundary value)
                dof_handler->boundary_dof_iterate(
                    MARS_LAMBDA(const ::mars::Integer local_dof) {
                        sp.vector_apply_constraints(local_dof, vec, value);
                    },
                    side);
            }

            void apply_zero_constraints(Vector &v, const std::string side) override {
                auto sp = *sparsity_pattern;
                auto vec = v.raw_type()->getLocalView<::mars::KokkosSpace>();
                // BC set values to constraint value to zero
                dof_handler->boundary_dof_iterate(
                    MARS_LAMBDA(const ::mars::Integer local_dof) { sp.apply_zero_constraints(local_dof, vec); }, side);
            }

            void copy_at_constrained_nodes(const Vector &in, Vector &out, const std::string side) override {
                auto sp = *sparsity_pattern;
                auto in_view = in.raw_type()->getLocalView<::mars::KokkosSpace>();
                auto out_view = out.raw_type()->getLocalView<::mars::KokkosSpace>();

                auto sp_dof_handler = sp.get_dof_handler();

                // BC set values to constraint value (i.e., boundary value)
                dof_handler->boundary_dof_iterate(
                    MARS_LAMBDA(const ::mars::Integer local_dof) {
                        // sp.vector_apply_constraints(local_dof, vec, value);
                        auto idx = sp_dof_handler.local_to_owned_index(local_dof);
                        out_view(idx, 0) = in_view(idx, 0);
                    },
                    side);
            }

            /* void system_apply_constraints(Matrix &m, Vector &v) override {
                vector_apply_constraints(v);
                matrix_apply_constraints(m, 1.0);
            } */

            SizeType n_local_dofs() override { return dof_handler->get_owned_dof_size(); };
            SizeType n_dofs() override { return dof_handler->get_global_dof_size(); };

            auto factory() -> Factory & override { return ConcreteFactory<DMesh>::instance(); };

            void init(DMesh &mesh_impl, int block_size) {
                dof_handler = std::make_shared<DofHandler>(&mesh_impl);  //, mesh->raw_type_context());
                dof_handler->set_block(block_size);
                dof_handler->enumerate_dofs();

                fe_dof_map = std::make_shared<FEDofMap>(build_fe_dof_map(*dof_handler));

                sparsity_pattern = std::make_shared<SPattern>(*dof_handler);
                sparsity_pattern->build_pattern(*fe_dof_map);
            }

            inline SPattern &get_sparsity_pattern() { return *sparsity_pattern; }
            inline DofHandler &get_dof_handler() { return *dof_handler; }
            inline FEDofMap &get_fe_dof_map() { return *fe_dof_map; }

            inline const SPattern &get_sparsity_pattern() const { return *sparsity_pattern; }
            inline const DofHandler &get_dof_handler() const { return *dof_handler; }
            inline const FEDofMap &get_fe_dof_map() const { return *fe_dof_map; }

        private:
            std::shared_ptr<SPattern> sparsity_pattern;
            std::shared_ptr<DofHandler> dof_handler;
            std::shared_ptr<FEDofMap> fe_dof_map;
        };

    }  // namespace mars
}  // namespace utopia
#endif

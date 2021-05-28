#ifndef UTOPIA_MARS_FEHANDLER_HPP
#define UTOPIA_MARS_FEHANDLER_HPP

#include "utopia_mars_FunctionSpace.hpp"
#include "mars.hpp"
#include "mars_base.hpp"
#include "mars_context.hpp"
#include "mars_distributed_dof_management.hpp"

#include "utopia_mars_Factory_impl.hpp"

namespace utopia {
    namespace mars {

        template <class DMesh, int Degree_ = 1>
        class FEHandler : public IFEHandler {
        public:
            static constexpr ::mars::Integer Degree = Degree_;
            using DofHandler = ::mars::DofHandler<DMesh, Degree>;
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

            void matrix_apply_constraints(Matrix &m, const Scalar diag_value) override {
                // BC set constrained rows to zero, except diagonal where you set diag_value
            }

            void vector_apply_constraints(Vector &v) override {
                // BC set values to constraint value (i.e., boundary value)
            }

            void apply_zero_constraints(Vector &vec) override {
                // BC set values to constraint value to zero
            }

            void system_apply_constraints(Matrix &m, Vector &v) override {
                vector_apply_constraints(v);
                matrix_apply_constraints(m, 1.0);
            }

            SizeType n_local_dofs() override { return dof_handler->get_owned_dof_size(); };
            SizeType n_dofs() override { return dof_handler->get_global_dof_size(); };

            auto factory() -> Factory & override { return ConcreteFactory<DMesh>::instance(); };

            void init(DMesh &mesh_impl) {
                dof_handler = std::make_shared<DofHandler>(&mesh_impl);  //, mesh->raw_type_context());
                dof_handler->enumerate_dofs();

                fe_dof_map = std::make_shared<FEDofMap>(build_fe_dof_map(*dof_handler));

                sparsity_pattern = std::make_shared<SPattern>(*dof_handler);
                sparsity_pattern->build_pattern(*fe_dof_map);
            }

        private:
            std::shared_ptr<SPattern> sparsity_pattern;
            std::shared_ptr<DofHandler> dof_handler;
            std::shared_ptr<FEDofMap> fe_dof_map;
        };

    }  // namespace mars
}  // namespace utopia
#endif

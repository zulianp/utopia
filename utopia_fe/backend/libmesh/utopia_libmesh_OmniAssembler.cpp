#include "utopia_libmesh_OmniAssembler.hpp"

#include "utopia_libmesh_FunctionSpace.hpp"

#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMaterial.hpp"

namespace utopia {
    namespace libmesh {

        class OmniAssembler::Impl {
        public:
            // New abstraction
            std::shared_ptr<libmesh::FunctionSpace> space;

            // Legacy abstractions
            using LegacyFunctionSpace = utopia::LibMeshFunctionSpace;
            using LegacyProductFunctionSpace = utopia::ProductFunctionSpace<LegacyFunctionSpace>;
            using LegacyMaterial = utopia::UIMaterial<LegacyProductFunctionSpace, Matrix, Vector>;
            using LegacyForcingFunction = utopia::UIForcingFunction<LegacyProductFunctionSpace, Vector>;
            using LegacyForcedMaterial = utopia::ForcedMaterial<Matrix, Vector>;

            std::shared_ptr<LegacyProductFunctionSpace> legacy_space;
            std::shared_ptr<Model<Matrix, Vector>> legacy_model;

            void update_legacy_mirror() {
                legacy_space = std::make_shared<LegacyProductFunctionSpace>();

                for (int s = 0; s < space->n_subspaces(); ++s) {
                    legacy_space->add_subspace(std::make_shared<LegacyFunctionSpace>(
                        make_ref(space->raw_type()), space->system_id(), (*space)[s].subspace_id()));
                }
            }
        };

        OmniAssembler::OmniAssembler(const std::shared_ptr<libmesh::FunctionSpace> &space)
            : impl_(utopia::make_unique<Impl>()) {
            impl_->space = space;
            impl_->update_legacy_mirror();
        }

        OmniAssembler::~OmniAssembler() = default;

        bool OmniAssembler::assemble(const Vector &x, Matrix &jacobian, Vector &fun) {
            if (!impl_->legacy_model) {
                return false;
            }

            return impl_->legacy_model->assemble_hessian_and_gradient(x, jacobian, fun);
        }

        void OmniAssembler::read(Input &in) {
            auto material = utopia::make_unique<Impl::LegacyMaterial>(*impl_->legacy_space);
            auto forcing_function = utopia::make_unique<Impl::LegacyForcingFunction>(*impl_->legacy_space);

            in.get("material", *material);
            in.get("forcing_functions", *forcing_function);

            assert(material->good());

            impl_->legacy_model =
                std::make_shared<Impl::LegacyForcedMaterial>(std::move(material), std::move(forcing_function));
        }

    }  // namespace libmesh
}  // namespace utopia

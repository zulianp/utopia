#include "utopia_libmesh_OmniAssembler.hpp"

#include "utopia_libmesh_FunctionSpace.hpp"

#include "utopia_UFlow.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMaterial.hpp"

namespace utopia {
    namespace libmesh {

        class OmniAssembler::Impl {
        public:
            using Scalar = Traits<libmesh::FunctionSpace>::Scalar;
            // New abstraction
            std::shared_ptr<libmesh::FunctionSpace> space;

            // Legacy abstractions
            using LegacyFunctionSpace = utopia::LibMeshFunctionSpace;
            using LegacyProductFunctionSpace = utopia::ProductFunctionSpace<LegacyFunctionSpace>;
            using LegacyMaterial = utopia::UIMaterial<LegacyProductFunctionSpace, Matrix, Vector>;
            using LegacyForcingFunction = utopia::UIForcingFunction<LegacyProductFunctionSpace, Vector>;
            using LegacyForcedMaterial = utopia::ForcedMaterial<Matrix, Vector>;
            using LegacyFlow = utopia::UFlow<LegacyFunctionSpace, Matrix, Vector>;

            std::shared_ptr<LegacyProductFunctionSpace> legacy_space;
            std::shared_ptr<Model<Matrix, Vector>> legacy_model;
            std::shared_ptr<LegacyForcingFunction> forcing_function;
            Scalar rescale{1.0};

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

            if (!impl_->legacy_model->assemble_hessian_and_gradient(x, jacobian, fun)) {
                return false;
            }

            Vector ff;
            if (!impl_->forcing_function->eval(x, ff)) {
                return false;
            }

            if (impl_->rescale != 1.0) {
                ff *= impl_->rescale;
            }

            fun -= ff;
            return true;
        }

        static bool is_flow(const std::string &name) {
            return (name.find("Flow") != std::string::npos) || (name.find("LaplaceOperator") != std::string::npos);
        }

        void OmniAssembler::read(Input &in) {
            std::string type;
            in.get("material", [&](Input &node) { node.get("type", type); });

            impl_->forcing_function = utopia::make_unique<Impl::LegacyForcingFunction>(*impl_->legacy_space);
            in.get("forcing_functions", *impl_->forcing_function);

            if (is_flow(type)) {
                // Flow problems
                auto material = utopia::make_unique<Impl::LegacyFlow>(impl_->legacy_space->subspace(0));
                in.get("material", *material);
                impl_->rescale = material->rescaling();

                impl_->legacy_model = std::move(material);
            } else {
                // Elasticity
                auto material = utopia::make_unique<Impl::LegacyMaterial>(*impl_->legacy_space);
                in.get("material", *material);
                assert(material->good());

                impl_->rescale = material->rescaling();
                impl_->legacy_model = std::move(material);
            }
        }

    }  // namespace libmesh
}  // namespace utopia

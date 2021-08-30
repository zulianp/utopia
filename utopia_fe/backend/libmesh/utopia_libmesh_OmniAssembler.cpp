#include "utopia_libmesh_OmniAssembler.hpp"

#include "utopia_libmesh_FunctionSpace.hpp"

#include "utopia_libmesh_Legacy.hpp"

#include "utopia_UFlow.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMaterial.hpp"

#include "utopia_libmesh_Transport.hpp"

namespace utopia {
    namespace libmesh {

        class AssemblerRegistry {
        public:
            using FEAssembler_t = utopia::FEAssembler<utopia::libmesh::FunctionSpace>;
            using FEAssemblerPtr_t = std::shared_ptr<FEAssembler_t>;

            FEAssemblerPtr_t find_assembler(const std::string &name) const {
                auto it = assemblers_.find(name);
                if (it == assemblers_.end()) {
                    // utopia::err() << "Unable to find assembler with name " << name << ".\n";
                    // assert(false);
                    return nullptr;
                } else {
                    return it->second();
                }
            }

            template <typename Assembler_t>
            void register_assembler(const std::string &name) {
                assemblers_[name] = []() -> FEAssemblerPtr_t { return std::make_shared<Assembler_t>(); };
            }

            AssemblerRegistry() { register_assemblers(); }

        private:
            std::map<std::string, std::function<FEAssemblerPtr_t()>> assemblers_;

            void register_assemblers() {
                register_assembler<Transport>("Transport");
                register_assembler<Mass>("Mass");
            }
        };

        class OmniAssembler::Impl {
        public:
            using Scalar = Traits<libmesh::FunctionSpace>::Scalar;
            // New abstraction
            std::shared_ptr<libmesh::FunctionSpace> space;
            std::shared_ptr<Environment<libmesh::FunctionSpace>> env;
            AssemblerRegistry registry;

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

            void update_legacy_mirror() { legacy_space = make_legacy(*space); }
        };

        void OmniAssembler::set_environment(const std::shared_ptr<Environment<libmesh::FunctionSpace>> &env) {
            impl_->env = env;
        }

        OmniAssembler::OmniAssembler(const std::shared_ptr<libmesh::FunctionSpace> &space)
            : impl_(utopia::make_unique<Impl>()) {
            impl_->space = space;
            impl_->update_legacy_mirror();
        }

        OmniAssembler::~OmniAssembler() = default;

        bool OmniAssembler::apply(const Vector &x, Vector &hessian_times_x) {
            assert(false);
            Utopia::Abort("IMPLEMENT ME");
            return false;
        }

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

        bool OmniAssembler::assemble(const Vector &x, Matrix &jacobian) {
            if (!impl_->legacy_model) {
                return false;
            }

            return impl_->legacy_model->assemble_hessian(x, jacobian);
        }

        bool OmniAssembler::assemble(const Vector &x, Vector &fun) {
            if (!impl_->legacy_model) {
                return false;
            }

            if (!impl_->legacy_model->assemble_gradient(x, fun)) {
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

        bool OmniAssembler::assemble(Matrix &jacobian) {
            if (!impl_->legacy_model) {
                return false;
            }

            return impl_->legacy_model->assemble_hessian(jacobian);
        }

        static bool is_flow(const std::string &name) {
            return (name.find("Flow") != std::string::npos) || (name.find("LaplaceOperator") != std::string::npos);
        }

        void OmniAssembler::read(Input &in) {
            std::string type;
            in.get("material", [&](Input &node) { node.get("type", type); });

            impl_->forcing_function = utopia::make_unique<Impl::LegacyForcingFunction>(*impl_->legacy_space);
            in.get("forcing_functions", *impl_->forcing_function);

            if (!type.empty()) {
                auto new_assembler = impl_->registry.find_assembler(type);
                if (new_assembler) {
                    new_assembler->set_space(impl_->space);
                    new_assembler->set_environment(impl_->env);

                    in.get("material", *new_assembler);
                    impl_->legacy_model = new_assembler;
                } else if (is_flow(type)) {
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
            } else if (impl_->legacy_model) {
                in.get("material", *impl_->legacy_model);
            }

            assert(impl_->legacy_model);
            if (!impl_->legacy_model) {
                Utopia::Abort();
            }
        }

        void OmniAssembler::set_time(const std::shared_ptr<SimulationTime> &) {
            utopia::err() << "[Warning] libmesh::OmniAssembler ingores time, for the moment!\n";
        }

        bool OmniAssembler::is_linear() const {
            if (impl_->legacy_model) {
                return impl_->legacy_model->is_linear();
            } else {
                assert(false);
                return true;
            }
        }

    }  // namespace libmesh
}  // namespace utopia

#include "utopia_mars_OmniAssembler.hpp"
#include "utopia_mars_Factory.hpp"

#include <memory>

namespace utopia {
    namespace mars {

        class OmniAssembler::Impl {
        public:
            using FEAssemblerPtr_t = std::shared_ptr<utopia::mars::FEAssembler>;

            std::shared_ptr<mars::FunctionSpace> space;
            std::vector<FEAssemblerPtr_t> assemblers;
        };

        OmniAssembler::OmniAssembler(const std::shared_ptr<mars::FunctionSpace> &space)
            : impl_(utopia::make_unique<Impl>()) {
            set_space(space);
        }

        OmniAssembler::~OmniAssembler() {}

        bool OmniAssembler::assemble(const Vector &x, Matrix &mat, Vector &vec) {
            UTOPIA_TRACE_REGION_BEGIN("OmniAssembler::assemble(x,mat,vec)");
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->assemble(x, mat, vec);
                if (!ok) {
                    break;
                }
            }

            UTOPIA_TRACE_REGION_END("OmniAssembler::assemble(x,mat,vec)");
            return ok;
        }

        bool OmniAssembler::assemble(const Vector &x, Matrix &mat) {
            UTOPIA_TRACE_REGION_BEGIN("OmniAssembler::assemble(x,mat)");
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->assemble(x, mat);
                if (!ok) {
                    break;
                }
            }

            UTOPIA_TRACE_REGION_END("OmniAssembler::assemble(x,mat)");
            return ok;
        }

        bool OmniAssembler::assemble(const Vector &x, Vector &vec) {
            UTOPIA_TRACE_REGION_BEGIN("OmniAssembler::assemble(x,vec)");
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->assemble(x, vec);
                if (!ok) {
                    break;
                }
            }

            UTOPIA_TRACE_REGION_END("OmniAssembler::assemble(x,vec)");
            return ok;
        }

        bool OmniAssembler::assemble(Matrix &mat) {
            UTOPIA_TRACE_REGION_BEGIN("OmniAssembler::assemble(mat)");
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->assemble(mat);
                if (!ok) {
                    break;
                }
            }

            UTOPIA_TRACE_REGION_END("OmniAssembler::assemble(mat)");
            return ok;
        }

        void OmniAssembler::set_space(const std::shared_ptr<FunctionSpace> &space) { impl_->space = space; }

        void OmniAssembler::read(Input &in) {
            in.get("material", [this](Input &node) {
                auto ass = impl_->space->factory().new_assembler(node);

                if (ass) {
                    ass->set_space(impl_->space);
                    ass->read(node);

                    impl_->assemblers.push_back(std::move(ass));
                }
            });

            if (impl_->assemblers.empty()) {
                Utopia::Abort("assemblers cannot be empty!");
            }
        }

        void OmniAssembler::set_environment(const std::shared_ptr<Environment> &env) {
            for (auto &a_ptr : impl_->assemblers) {
                a_ptr->set_environment(env);
            }
        }

        bool OmniAssembler::is_linear() const {
            for (auto &a_ptr : impl_->assemblers) {
                if (!a_ptr->is_linear()) {
                    return false;
                }
            }

            return true;
        }

        std::string OmniAssembler::name() const { return "OmniAssembler"; }

        bool OmniAssembler::apply(const Vector &x, Vector &hessian_times_x) {
            if (!is_operator()) return false;

            UTOPIA_TRACE_REGION_BEGIN("OmniAssembler::apply(x,hessian_times_x)");
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->apply(x, hessian_times_x);
                if (!ok) {
                    break;
                }
            }

            UTOPIA_TRACE_REGION_END("OmniAssembler::apply(x,hessian_times_x)");
            return ok;
        }

        bool OmniAssembler::is_operator() const {
            for (auto &a_ptr : impl_->assemblers) {
                if (!a_ptr->is_operator()) {
                    return false;
                }
            }

            return true;
        }

        bool OmniAssembler::assemble(Vector &fun) {
            Utopia::Abort("IMPLEMENT ME");
            return false;
        }

        std::shared_ptr<OmniAssembler::Environment> OmniAssembler::environment() const {
            Utopia::Abort("IMPLEMENT ME");
            return nullptr;
        }

        std::shared_ptr<FunctionSpace> OmniAssembler::space() const {
            Utopia::Abort("IMPLEMENT ME");
            return nullptr;
        }

        void OmniAssembler::set_time(const std::shared_ptr<SimulationTime> &time) { Utopia::Abort("IMPLEMENT ME"); }

    }  // namespace mars
}  // namespace utopia

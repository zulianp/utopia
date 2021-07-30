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
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->assemble(x, mat, vec);
                if (!ok) {
                    return false;
                }
            }

            return ok;
        }

        bool OmniAssembler::assemble(const Vector &x, Matrix &mat) {
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->assemble(x, mat);
                if (!ok) {
                    return false;
                }
            }

            return ok;
        }

        bool OmniAssembler::assemble(const Vector &x, Vector &vec) {
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->assemble(x, vec);
                if (!ok) {
                    return false;
                }
            }

            return ok;
        }

        bool OmniAssembler::assemble(Matrix &mat) {
            bool ok = false;
            for (auto &a_ptr : impl_->assemblers) {
                ok = a_ptr->assemble(mat);
                if (!ok) {
                    return false;
                }
            }

            return ok;
        }

        void OmniAssembler::set_space(const std::shared_ptr<FunctionSpace> &space) { impl_->space = space; }

        void OmniAssembler::read(Input &in) {
            auto ass = impl_->space->factory().new_assembler(in);

            if (ass) {
                ass->set_space(impl_->space);
                ass->read(in);

                impl_->assemblers.push_back(std::move(ass));
            }

            if (impl_->assemblers.empty()) {
                Utopia::Abort("assemblers cannot be empty!");
            }
        }

        void OmniAssembler::set_environment(const std::shared_ptr<Environment<mars::FunctionSpace>> &env) {
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

    }  // namespace mars
}  // namespace utopia

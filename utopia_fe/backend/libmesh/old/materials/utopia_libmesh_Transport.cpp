#include "utopia_libmesh_Transport.hpp"

#include "utopia_make_unique.hpp"

#include "utopia_libmesh_Legacy.hpp"

/// Legacy includs
#include "utopia_Flow.hpp"
#include "utopia_StabilizeTransport.hpp"
#include "utopia_UIScalarSampler.hpp"

namespace utopia {
    namespace libmesh {

        class Transport::Impl {
        public:
            inline bool initialized() const { return static_cast<bool>(legacy_space); }

            std::shared_ptr<LegacyProductFunctionSpace> legacy_space;
            std::shared_ptr<Field> pressure;
            UIScalarFunction<Scalar> diffusion_function;
            Scalar coeff{1.0};
            bool stabilize_transport{true};
            bool verbose{false};
        };

        Transport::Transport() : impl_(utopia::make_unique<Impl>()) {}
        Transport::~Transport() = default;

        void Transport::init() {
            auto s = this->space();
            assert(s);

            if (s) {
                impl_->legacy_space = make_legacy(*s);
            }
        }

        void Transport::set_pressure_field(const std::shared_ptr<Field> &field) { impl_->pressure = field; }

        bool Transport::valid() const { return impl_->legacy_space && impl_->pressure; }

        bool Transport::assemble(const Vector &x, Matrix &jacobian, Vector &fun) {
            if (!impl_->initialized()) {
                init();
            }

            assert(valid());
            if (!valid()) return false;

            auto &pressure_w = impl_->pressure->data();
            auto &C = impl_->legacy_space->subspace(0);
            auto &mesh = C.mesh();
            auto &dof_map = C.dof_map();
            const int dim = mesh.spatial_dimension();

            auto c = trial(C);
            auto q = test(C);
            auto ph = interpolate(pressure_w, c);

            auto sampler_fun = ctx_fun(impl_->diffusion_function.sampler());
            auto vel = sampler_fun * impl_->coeff * grad(ph);
            auto b_form = (inner(inner(-grad(c), vel), q) * dX);

            utopia::assemble(b_form, jacobian);

            if (impl_->verbose) {
                Scalar norm1_jac = norm1(jacobian);
                utopia::out() << "norm1 Jacobian (no stab): " << norm1_jac << '\n';
            }

            if (impl_->stabilize_transport) {
                Matrix jacobian_out;
                transport_stabilization(jacobian, jacobian_out);
                rename(jacobian.name(), jacobian_out);
                jacobian = std::move(jacobian_out);
            }

            if (impl_->verbose) {
                Scalar norm1_jac = norm1(jacobian);
                utopia::out() << "norm1 Jacobian (with stab): " << norm1_jac << '\n';
            }

            return true;
        }

        void Transport::clear() {}

        void Transport::read(Input &in) {
            Super::read(in);
            auto env = this->environment();

            if (env) {
                std::string pressure_field;
                in.get("pressure_field", pressure_field);
                auto p = env->find_field(*this->space(), pressure_field);

                assert(p);

                if (p) {
                    impl_->pressure = p;
                }
            } else {
                std::string pressure_field;
                in.get("pressure_field", pressure_field);

                if (!pressure_field.empty()) {
                    assert(false);
                    Utopia::Abort("In order to retrive the pressure_field, The env must be defined!");
                }
            }

            in.get("diffusion_function", impl_->diffusion_function);
            in.get("coeff", impl_->coeff);
            in.get("stabilize_transport", impl_->stabilize_transport);
            in.get("verbose", impl_->verbose);

            init();

            if (impl_->verbose) {
                utopia::out() << "-----------------------------\n";
                utopia::out() << "Transport\n";
                if (impl_->pressure) {
                    utopia::out() << "Pressure field:\t" << impl_->pressure->name() << '\n';
                }

                utopia::out() << "coeff:\t" << impl_->coeff << '\n';
                utopia::out() << "stabilize_transport:\t" << impl_->stabilize_transport << '\n';
                utopia::out() << "diffusion_function: ";
                impl_->diffusion_function.describe(utopia::out().stream());
                utopia::out() << "-----------------------------\n";
            }
        }

        /////////////////////////////////////////

        class Mass::Impl {
        public:
            inline bool initialized() const { return static_cast<bool>(legacy_space); }

            std::shared_ptr<LegacyProductFunctionSpace> legacy_space;
            // std::shared_ptr<Field> pressure;
            UIScalarFunction<Scalar> density_function;
            Scalar density{1.0};
            bool lumped{true};
            bool verbose{false};
        };

        Mass::Mass() : impl_(utopia::make_unique<Impl>()) {}
        Mass::~Mass() = default;

        void Mass::init() {
            auto s = this->space();
            assert(s);

            if (s) {
                impl_->legacy_space = make_legacy(*s);
            }
        }

        bool Mass::valid() const { return impl_->initialized(); }

        bool Mass::assemble(const Vector &x, Matrix &jacobian, Vector &fun) {
            if (!impl_->initialized()) {
                init();
            }

            assert(valid());
            if (!valid()) return false;

            auto &space = *impl_->legacy_space;
            auto &mesh = space[0].mesh();
            auto &dof_map = space[0].dof_map();
            const int dim = mesh.spatial_dimension();

            auto b_form =
                impl_->density * inner(ctx_fun(impl_->density_function.sampler()) * trial(space), test(space)) * dX;

            utopia::assemble(b_form, jacobian);

            if (impl_->lumped) {
                Vector mass_vector = sum(jacobian, 1);
                jacobian = diag(mass_vector);
            } else {
                utopia::assemble(b_form, jacobian);
            }

            fun = jacobian * x;
            return true;
        }

        void Mass::clear() {}

        void Mass::read(Input &in) {
            Super::read(in);
            in.get("density_function", impl_->density_function);
            in.get("density", impl_->density);
            in.get("lumped", impl_->lumped);
            in.get("verbose", impl_->verbose);

            init();

            if (impl_->verbose) {
                utopia::out() << "-----------------------------\n";
                utopia::out() << "Mass\n";
                utopia::out() << "lumped:\t" << impl_->lumped << '\n';
                utopia::out() << "density:\t" << impl_->density << '\n';
                utopia::out() << "density_function: ";
                impl_->density_function.describe(utopia::out().stream());
                utopia::out() << "-----------------------------\n";
            }
        }

    }  // namespace libmesh
}  // namespace utopia

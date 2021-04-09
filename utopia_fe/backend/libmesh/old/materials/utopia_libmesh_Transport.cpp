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
            UIScalarFunction<Scalar> diffusion_sampler;
            Scalar coeff{1.0};
            bool stabilize_transport{true};
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

            auto sampler_fun = ctx_fun(impl_->diffusion_sampler.sampler());
            auto vel = sampler_fun * impl_->coeff * grad(ph);
            auto b_form = (inner(inner(-grad(c), vel), q) * dX);

            utopia::assemble(b_form, jacobian);

            if (impl_->stabilize_transport) {
                Matrix jacobian_out;
                transport_stabilization(jacobian, jacobian_out);
                jacobian = std::move(jacobian_out);
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
                auto p = env->find_field(pressure_field);

                assert(p);

                if (p) {
                    impl_->pressure = p;
                }
            }

            in.get("diffusion_sampler", impl_->diffusion_sampler);
            in.get("coeff", impl_->coeff);
            in.get("stabilize_transport", impl_->stabilize_transport);

            init();
        }

    }  // namespace libmesh
}  // namespace utopia

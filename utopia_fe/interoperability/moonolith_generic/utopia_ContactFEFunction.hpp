#ifndef UTOPIA_CONTACTF_EF_UNCTION_HPP
#define UTOPIA_CONTACTF_EF_UNCTION_HPP

#include "utopia_BoxConstrainedFEFunction.hpp"
// #include "utopia_IContact.hpp"
#include "utopia_IsNotSupported.hpp"
#include "utopia_fe_Core.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ContactFEFunction final : public BoxConstrainedFEFunction<FunctionSpace> {
    public:
        using Super = utopia::BoxConstrainedFEFunction<FunctionSpace>;
        using Vector_t = typename Traits<FunctionSpace>::Vector;
        using Matrix_t = typename Traits<FunctionSpace>::Matrix;

        using Size_t = typename Traits<FunctionSpace>::SizeType;
        using Scalar_t = typename Traits<FunctionSpace>::Scalar;
        using Communicator_t = typename Traits<FunctionSpace>::Communicator;
        using Mesh_t = typename Traits<FunctionSpace>::Mesh;
        using Contact_t = utopia::Contact<FunctionSpace>;

        bool has_nonlinear_constraints() const override { return !linear_contact_; }

        inline std::shared_ptr<Vector_t> selection() override {
            auto s = std::make_shared<Vector_t>(contact_->is_contact());
            return s;
        }

        bool update_constraints(const Vector_t &x) {
            this->space()->displace(x);
            bool ok = contact_->assemble(*this->space());

            if (debug_) {
                static int iter_debug = 0;
                ouput_debug_data(iter_debug++, x);
            }

            this->space()->displace(-x);
            return ok;
        }

        bool has_transformation() const override { return true; }

        void transform(const Vector_t &x, Vector_t &x_constrained) override { contact_->transform(x, x_constrained); }

        void inverse_transform(const Vector_t &x_constrained, Vector_t &x) override {
            contact_->inverse_transform(x_constrained, x);
        }

        void transform(const Matrix_t &H, Matrix_t &H_constrained) override { contact_->transform(H, H_constrained); }

        bool constraints_gradient(const Vector_t &x, BoxConstraints<Vector_t> &box) override {
            update_constraints(x);

            if (!box.upper_bound()) {
                box.upper_bound() = std::make_shared<Vector_t>();
            }

            (*box.upper_bound()) = contact_->gap();
            return true;
        }

        ContactFEFunction(const std::shared_ptr<FEFunctionInterface<FunctionSpace>> &unconstrained)
            : Super(unconstrained) {}

        void ouput_debug_data(const Size_t iter_debug, const Vector_t &) const {
            // Output extras
            {
                Vector_t gap_zeroed = e_mul(contact_->is_contact(), contact_->gap());
                Vector_t rays = contact_->normals();

                {
                    auto rays_view = local_view_device(rays);
                    auto gap_view = const_local_view_device(contact_->gap());

                    int dim = this->space()->mesh().spatial_dimension();

                    parallel_for(
                        local_range_device(rays), UTOPIA_LAMBDA(const Size_t i) {
                            Size_t block = i / dim;
                            Size_t g_idx = dim * block;

                            auto ri = rays_view.get(i);
                            ri *= gap_view.get(g_idx);
                            rays_view.set(i, ri);
                        });
                }

                this->space()->write("rays_" + std::to_string(iter_debug) + ".e", rays);
            }
        }

        void read(Input &in) override {
            Super::read(in);

            in.get("linear_contact", linear_contact_);
            in.get("debug", debug_);

            if (!contact_) {
                contact_ = std::make_shared<Contact_t>();
                in.get("contact", [&](Input &node) { contact_->read(node); });
            }
        }

    private:
        std::shared_ptr<Contact_t> contact_;
        bool linear_contact_{false};
        bool debug_{false};
    };

}  // namespace utopia

#endif  // UTOPIA_CONTACT_FE_FUNCTION_HPP

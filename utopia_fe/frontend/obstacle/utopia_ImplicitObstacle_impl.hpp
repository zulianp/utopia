#ifndef UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP
#define UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP

#include "utopia_IPTransfer.hpp"
#include "utopia_ImplicitObstacle.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ImplicitObstacle<FunctionSpace>::Impl {
    public:
        std::shared_ptr<FunctionSpace> domain;
        std::shared_ptr<Transfer> transfer;
        std::shared_ptr<Field> domain_distance;
        std::shared_ptr<Field> gap;
        std::shared_ptr<GradientField> normals;
        Matrix orthogonal_trafo;

        bool update_transfer{true};
    };

    template <class FunctionSpace>
    ImplicitObstacle<FunctionSpace>::ImplicitObstacle() : impl_(utopia::make_unique<Impl>()) {}

    template <class FunctionSpace>
    ImplicitObstacle<FunctionSpace>::~ImplicitObstacle() = default;

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::read(Input &in) {
        if (!impl_->domain) {
            impl_->domain = utopia::make_unique<FunctionSpace>();
            impl_->domain_distance = utopia::make_unique<Field>();
            in.require("domain", [this](Input &node) {
                // Must have a field associated with the space
                impl_->domain->read_with_state(node, *impl_->domain_distance);
            });
        }

        in.get("update_transfer", impl_->update_transfer);
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::describe(std::ostream &os) const {}

    template <class FunctionSpace>
    bool ImplicitObstacle<FunctionSpace>::init(const std::shared_ptr<FunctionSpace> &domain) {
        impl_->domain = domain;
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::set_transfer(const std::shared_ptr<Transfer> &transfer) {
        impl_->transfer = transfer;
    }

    template <class FunctionSpace>
    bool ImplicitObstacle<FunctionSpace>::assemble(FunctionSpace &space) {
        if (!impl_->transfer || impl_->update_transfer) {
            FETransfer<FunctionSpace> transfer;
            transfer.init(impl_->domain, make_ref(space));
            impl_->transfer = transfer.template build_transfer<IPTransfer<Matrix, Vector>>();
        }

        auto gap = std::make_shared<Field>();
        auto normals = std::make_shared<GradientField>();

        space.create_field(*gap);
        normals->init(*gap);
        normals->normalize();

        impl_->gap = gap;
        impl_->normals = normals;
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::transform(const Matrix &in, Matrix &out) {}

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::transform(const Vector &in, Vector &out) {}

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::inverse_transform(const Vector &in, Vector &out) {}

    template <class FunctionSpace>
    const typename ImplicitObstacle<FunctionSpace>::Vector &ImplicitObstacle<FunctionSpace>::gap() const {
        assert(impl_->gap);
        return impl_->gap->data();
    }

    template <class FunctionSpace>
    const typename ImplicitObstacle<FunctionSpace>::Vector &ImplicitObstacle<FunctionSpace>::is_contact() const {}

    template <class FunctionSpace>
    const typename ImplicitObstacle<FunctionSpace>::Vector &ImplicitObstacle<FunctionSpace>::normals() const {
        assert(impl_->normals);
        return impl_->normals->data();
    }

}  // namespace utopia

#endif  // UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP

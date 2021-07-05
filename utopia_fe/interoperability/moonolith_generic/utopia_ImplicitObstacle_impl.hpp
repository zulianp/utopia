#ifndef UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP
#define UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP

#include "utopia.hpp"
#include "utopia_IPTransfer.hpp"

#include "utopia_Field.hpp"
#include "utopia_ImplicitObstacle.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_moonolith_HouseholderReflection.hpp"

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
        Vector is_contact;

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

                Scalar norm_dd = norm2(impl_->domain_distance->data());
                utopia::out() << "norm_dd: " << norm_dd << "\n";
            });
        }

        in.get("update_transfer", impl_->update_transfer);
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::describe(std::ostream &os) const {
        UTOPIA_UNUSED(os);
    }

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

        impl_->transfer->interpolate(impl_->domain_distance->data(), gap->data());
        normals->init_and_normalize(*gap);

        impl_->gap = gap;
        impl_->normals = normals;

        impl_->is_contact.zeros(layout(normals->data()));

        int n_var = space.n_var();

        auto r = local_range_device(impl_->is_contact);
        RangeDevice<Vector> rd(r.begin(), r.begin() + r.extent() / n_var);
        auto view = local_view_device(impl_->is_contact);

        parallel_for(
            rd, UTOPIA_LAMBDA(const SizeType i) { view.set(i * n_var, 1.); });

        // Comment me out
        {
            utopia::out() << "n_dofs: " << space.n_dofs() << " n_nodes: " << space.mesh().n_nodes()
                          << " nn: " << impl_->normals->data().size() << "\n";

            Field normals_out("normals");
            space.create_nodal_vector_field(space.mesh().spatial_dimension(), normals_out);
            normals_out.data() = impl_->normals->data();
            space.backend_set_nodal_field(normals_out);

            IO<FunctionSpace> io(space);
            io.set_output_path("implicit_obstacle_out.e");
            io.register_output_field(normals_out.name());
            io.write(this->gap(), 0, 0);
        }

        int spatial_dimension = space.mesh().spatial_dimension();

        switch (spatial_dimension) {
            case 2: {
                HouseholderReflectionForContact<Matrix, Vector, 2>::build(
                    impl_->is_contact, normals->data(), impl_->orthogonal_trafo);
                break;
            }
            case 3: {
                HouseholderReflectionForContact<Matrix, Vector, 3>::build(
                    impl_->is_contact, normals->data(), impl_->orthogonal_trafo);
                break;
            }
            default: {
                Utopia::Abort("Unsupported dimension!");
                return false;
            }
        }

        return true;
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::transform(const Matrix &in, Matrix &out) {
        out = impl_->orthogonal_trafo * in;
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::transform(const Vector &in, Vector &out) {
        out = impl_->orthogonal_trafo * in;
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::inverse_transform(const Vector &in, Vector &out) {
        out = impl_->orthogonal_trafo * in;
    }

    template <class FunctionSpace>
    const typename ImplicitObstacle<FunctionSpace>::Vector &ImplicitObstacle<FunctionSpace>::gap() const {
        assert(impl_->gap);
        return impl_->gap->data();
    }

    template <class FunctionSpace>
    const typename ImplicitObstacle<FunctionSpace>::Vector &ImplicitObstacle<FunctionSpace>::is_contact() const {
        return impl_->is_contact;
    }

    template <class FunctionSpace>
    const typename ImplicitObstacle<FunctionSpace>::Vector &ImplicitObstacle<FunctionSpace>::normals() const {
        assert(impl_->normals);
        return impl_->normals->data();
    }

}  // namespace utopia

#endif  // UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP

#ifndef UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP
#define UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP

#include "utopia.hpp"
#include "utopia_FETransferOptions.hpp"
#include "utopia_Field.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_ImplicitObstacle.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_moonolith_HouseholderReflection.hpp"

namespace utopia {

    template <class FunctionSpace>
    class ImplicitObstacle<FunctionSpace>::Impl {
    public:
        // Domain quantities
        std::shared_ptr<FunctionSpace> domain;
        std::shared_ptr<Field> domain_distance;
        std::shared_ptr<GradientField> domain_gradients;

        std::shared_ptr<Transfer> transfer;
        std::shared_ptr<Field> gap;
        std::shared_ptr<GradientField> normals;

        Matrix orthogonal_trafo;
        Vector is_contact;

        Scalar infinity{100000};

        bool update_transfer{true};
        bool export_tensors{false};
        bool shift_field{false};
        bool volume_to_surface{false};
        bool has_covering{false};
        Scalar field_rescale{1.0};
        Scalar field_offset{0.0};
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

        in.get("shift_field", impl_->shift_field);
        in.get("update_transfer", impl_->update_transfer);
        in.get("export_tensors", impl_->export_tensors);
        in.get("infinity", impl_->infinity);
        in.get("field_rescale", impl_->field_rescale);
        in.get("field_offset", impl_->field_offset);
        in.get("volume_to_surface", impl_->volume_to_surface);
        in.get("has_covering", impl_->has_covering);

        if (impl_->shift_field) {
            Scalar min_dd = min(impl_->domain_distance->data());
            Scalar offset = impl_->field_offset;

            utopia::out() << "min_dd: " << min_dd << "\n";

            impl_->domain_distance->data().transform_values(
                UTOPIA_LAMBDA(const Scalar &val)->Scalar { return val - min_dd + offset; });

            if (impl_->field_rescale != 1.0) {
                impl_->domain_distance->data() *= impl_->field_rescale;
            }
        }

        compute_gradients();
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::compute_gradients() {
        impl_->domain_gradients = std::make_shared<GradientField>();
        impl_->domain_gradients->init_and_normalize(*impl_->domain_distance);
        impl_->domain_gradients->data() *= -1;

        // Comment me out
        {
            utopia::out() << "n_var: " << impl_->domain->n_var() << " "
                          << "n_dofs: " << impl_->domain->n_dofs() << " n_nodes: " << impl_->domain->mesh().n_nodes()
                          << " tensor_size: " << impl_->domain_gradients->tensor_size() << " "
                          << " nn: " << impl_->domain_gradients->data().size() << "\n";

            Field normals_out("normals");
            impl_->domain->create_nodal_vector_field(impl_->domain->mesh().spatial_dimension(), normals_out);
            normals_out.data() = impl_->domain_gradients->data();
            normals_out.set_tensor_size(impl_->domain->mesh().spatial_dimension());
            impl_->domain->backend_set_nodal_field(normals_out);

            IO<FunctionSpace> io(*impl_->domain);
            io.set_output_path("implicit_obstacle_out.e");
            io.register_output_field(normals_out.name());
            io.write(impl_->domain_distance->data(), 0, 0);
        }
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::describe(std::ostream &os) const {
        UTOPIA_UNUSED(os);
    }

    template <class FunctionSpace>
    bool ImplicitObstacle<FunctionSpace>::init(const std::shared_ptr<FunctionSpace> &domain) {
        impl_->domain = domain;
        return true;
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::set_transfer(const std::shared_ptr<Transfer> &transfer) {
        impl_->transfer = transfer;
    }

    template <class FunctionSpace>
    bool ImplicitObstacle<FunctionSpace>::assemble(FunctionSpace &space) {
        FETransferOptions opts;
        opts.has_covering = false;

        FETransfer<FunctionSpace> transfer;
        transfer.set_options(opts);

        InputParameters params;
        if (impl_->volume_to_surface) {
            params.set("volume_to_surface", true);
        }

        params.set("has_covering", impl_->has_covering);
        transfer.read(params);

        transfer.init(impl_->domain, make_ref(space));
        impl_->transfer = transfer.template build_transfer<IPTransfer<Matrix, Vector>>();

        auto gap = std::make_shared<Field>();
        auto normals = std::make_shared<GradientField>();
        int dim = space.mesh().spatial_dimension();

        space.create_field(*gap);
        space.create_nodal_vector_field(dim, *normals);

        transfer.apply(impl_->domain_distance->data(), gap->data());
        transfer.apply(impl_->domain_gradients->data(), normals->data());

        normals->normalize();
        assert(!normals->data().has_nan_or_inf());

        impl_->gap = gap;
        impl_->normals = normals;
        impl_->is_contact.zeros(layout(normals->data()));

        int n_var = space.n_var();

        auto r = local_range_device(impl_->is_contact);
        RangeDevice<Vector> rd(r.begin(), r.begin() + r.extent() / n_var);

        auto gap_view = local_view_device(gap->data());
        auto normals_view = local_view_device(normals->data());
        Scalar infty = impl_->infinity;

        {
            auto is_contact_view = local_view_device(impl_->is_contact);

            parallel_for(
                rd, UTOPIA_LAMBDA(const SizeType i) {
                    Scalar norm_n = 0.0;
                    for (int d = 0; d < n_var; ++d) {
                        auto x = normals_view.get(i * n_var + d);
                        norm_n += x * x;
                    }

                    int start = 0;
                    if (norm_n != 0.) {
                        start = 1;
                        is_contact_view.set(i * n_var, 1.);
                    }

                    for (int k = start; k < n_var; ++k) {
                        gap_view.set(i * n_var + k, infty);
                    }
                });
        }

        Vector ones(layout(impl_->is_contact), 1);
        space.apply_zero_constraints(ones);

        {
            auto one_view = local_view_device(ones);
            auto is_contact_view = local_view_device(impl_->is_contact);

            // Remove Dirichlet
            parallel_for(
                rd, UTOPIA_LAMBDA(const SizeType i) {
                    bool is_dirichlet = false;
                    for (int k = 0; k < n_var; ++k) {
                        if (one_view.get(i * n_var + k) < 0.99) {
                            is_dirichlet = true;
                        }
                    }

                    bool is_c = (!is_dirichlet) && (is_contact_view.get(i * n_var) > 0.99);

                    if (!is_c) {
                        for (int k = 0; k < n_var; ++k) {
                            gap_view.set(i * n_var + k, infty);
                            is_contact_view.set(i * n_var + k, 0);
                        }
                    }
                });
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

        if (impl_->export_tensors) {
            rename("n", impl_->normals->data());
            rename("O", impl_->orthogonal_trafo);
            rename("ind", impl_->is_contact);
            rename("g", impl_->gap->data());

            write("load_n.m", this->normals());
            write("load_O.m", impl_->orthogonal_trafo);
            write("load_ind.m", this->is_contact());
            write("load_g.m", this->gap());
        }

        return true;
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::transform(const Matrix &in, Matrix &out) {
        out = transpose(impl_->orthogonal_trafo) * in * impl_->orthogonal_trafo;
    }

    template <class FunctionSpace>
    void ImplicitObstacle<FunctionSpace>::transform(const Vector &in, Vector &out) {
        out = transpose(impl_->orthogonal_trafo) * in;
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

    template <class FunctionSpace>
    std::shared_ptr<typename Traits<FunctionSpace>::Matrix> ImplicitObstacle<FunctionSpace>::orthogonal_transformation()
        const {
        return make_ref(impl_->orthogonal_trafo);
    }

}  // namespace utopia

#endif  // UTOPIA_IMPLICIT_OBSTACLE_IMPL_HPP

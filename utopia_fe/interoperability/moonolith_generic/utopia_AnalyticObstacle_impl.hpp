#ifndef UTOPIA_ANALITYC_OBSTACLE_IMPL_HPP
#define UTOPIA_ANALITYC_OBSTACLE_IMPL_HPP

#include "utopia.hpp"
#include "utopia_IPTransfer.hpp"
#include "utopia_SymbolicFunction.hpp"

#include "utopia_AnalyticObstacle.hpp"
#include "utopia_Field.hpp"
#include "utopia_fe_Core.hpp"
#include "utopia_moonolith_HouseholderReflection.hpp"

namespace utopia {

    template <class FunctionSpace>
    class AnalyticObstacle<FunctionSpace>::Impl {
    public:
        std::shared_ptr<Field> gap;
        std::shared_ptr<GradientField> normals;
        Matrix orthogonal_trafo;
        Vector is_contact;
        Scalar infinity{1e10};
        bool export_tensors{false};

        class Fun : public Configurable {
        public:
            utopia::SymbolicFunction x;
            utopia::SymbolicFunction y;
            utopia::SymbolicFunction z;
            utopia::SymbolicFunction normal_x;
            utopia::SymbolicFunction normal_y;
            utopia::SymbolicFunction normal_z;

            std::string boundary_name;

            Fun() : x("0"), y("0"), z("0"), normal_x("0"), normal_y("-1"), normal_z("0") {}

            void read(Input &in) override {
                std::string expr_x = "0";
                std::string expr_y = "0";
                std::string expr_z = "0";

                std::string expr_normal_x = "0";
                std::string expr_normal_y = "-1";
                std::string expr_normal_z = "0";

                in.get("x", expr_x);
                in.get("y", expr_y);
                in.get("z", expr_z);
                in.get("normal_x", expr_normal_x);
                in.get("normal_y", expr_normal_y);
                in.get("normal_z", expr_normal_z);

                x = utopia::SymbolicFunction(expr_x);
                y = utopia::SymbolicFunction(expr_y);
                z = utopia::SymbolicFunction(expr_z);

                normal_x = utopia::SymbolicFunction(expr_normal_x);
                normal_y = utopia::SymbolicFunction(expr_normal_y);
                normal_z = utopia::SymbolicFunction(expr_normal_z);

                in.get("boundary_name", boundary_name);
            }
        };

        Fun fun;

        bool update_transfer{true};
    };

    template <class FunctionSpace>
    AnalyticObstacle<FunctionSpace>::AnalyticObstacle() : impl_(utopia::make_unique<Impl>()) {}

    template <class FunctionSpace>
    AnalyticObstacle<FunctionSpace>::~AnalyticObstacle() = default;

    template <class FunctionSpace>
    void AnalyticObstacle<FunctionSpace>::read(Input &in) {
        in.get("fun", impl_->fun);
        in.get("export_tensors", impl_->export_tensors);
    }

    template <class FunctionSpace>
    void AnalyticObstacle<FunctionSpace>::describe(std::ostream &os) const {
        UTOPIA_UNUSED(os);
    }

    template <class FunctionSpace>
    bool AnalyticObstacle<FunctionSpace>::init(const std::shared_ptr<FunctionSpace> &domain) {
        impl_->domain = domain;
    }

    template <class FunctionSpace>
    void AnalyticObstacle<FunctionSpace>::set_transfer(const std::shared_ptr<Transfer> &transfer) {
        impl_->transfer = transfer;
    }

    template <class FunctionSpace>
    bool AnalyticObstacle<FunctionSpace>::assemble(FunctionSpace &space) {
        auto gap = std::make_shared<Field>();
        auto normals = std::make_shared<GradientField>();
        int dim = space.mesh().spatial_dimension();

        space.create_field(*gap);
        space.create_field(*normals);

        impl_->gap = gap;
        impl_->normals = normals;
        impl_->is_contact.zeros(layout(normals->data()));

        impl_->gap->data().set(impl_->infinity);

        int n_var = space.n_var();

        {
            // Write<Vector> wg(impl_->gap->data()), wn(impl_->normals->data());

            auto g_view = view_device(impl_->gap->data());
            auto n_view = view_device(impl_->normals->data());
            auto c_view = local_view_device(impl_->is_contact);

            // FIXME (will never work on GPU)
            space.node_eval([this, dim, g_view, n_view, c_view, n_var](const SizeType idx, const Scalar *point) {
                Scalar coords[4] = {0.0, 0.0, 0.0, 0.0};
                Scalar normal[3] = {0.0, 0.0, 0.0};
                Scalar p[3] = {0.0, 0.0, 0.0};

                for (int d = 0; d < dim; ++d) {
                    coords[d] = point[d];
                }

                p[0] = impl_->fun.x.eval(coords[0], coords[1], coords[2], coords[3]);
                p[1] = impl_->fun.y.eval(coords[0], coords[1], coords[2], coords[3]);
                p[2] = impl_->fun.z.eval(coords[0], coords[1], coords[2], coords[3]);

                normal[0] = impl_->fun.normal_x.eval(coords[0], coords[1], coords[2], coords[3]);
                normal[1] = impl_->fun.normal_y.eval(coords[0], coords[1], coords[2], coords[3]);
                normal[2] = impl_->fun.normal_z.eval(coords[0], coords[1], coords[2], coords[3]);

                Scalar g = 0.0;

                for (int d = 0; d < dim; ++d) {
                    auto x = (coords[d] - p[d]) * normal[d];
                    g += x;
                }

                g_view.set(idx * n_var, g);
                c_view.set(idx * n_var, 1);

                for (int d = 0; d < dim; ++d) {
                    // Invert the normal
                    n_view.set(idx * n_var + d, -normal[d]);
                }
            });
        }

        space.apply_zero_constraints(impl_->is_contact);
        space.apply_constraints(impl_->gap->data());
        space.apply_zero_constraints(normals->data());

        switch (dim) {
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

        // space.write("is_contact.e", impl_->is_contact);

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
    void AnalyticObstacle<FunctionSpace>::transform(const Matrix &in, Matrix &out) {
        out = impl_->orthogonal_trafo * in;
    }

    template <class FunctionSpace>
    void AnalyticObstacle<FunctionSpace>::transform(const Vector &in, Vector &out) {
        out = transpose(impl_->orthogonal_trafo) * in;
    }

    template <class FunctionSpace>
    void AnalyticObstacle<FunctionSpace>::inverse_transform(const Vector &in, Vector &out) {
        out = impl_->orthogonal_trafo * in;
    }

    template <class FunctionSpace>
    const typename AnalyticObstacle<FunctionSpace>::Vector &AnalyticObstacle<FunctionSpace>::gap() const {
        assert(impl_->gap);
        return impl_->gap->data();
    }

    template <class FunctionSpace>
    const typename AnalyticObstacle<FunctionSpace>::Vector &AnalyticObstacle<FunctionSpace>::is_contact() const {
        return impl_->is_contact;
    }

    template <class FunctionSpace>
    const typename AnalyticObstacle<FunctionSpace>::Vector &AnalyticObstacle<FunctionSpace>::normals() const {
        assert(impl_->normals);
        return impl_->normals->data();
    }

}  // namespace utopia

#endif  // UTOPIA_ANALITYC_OBSTACLE_IMPL_HPP

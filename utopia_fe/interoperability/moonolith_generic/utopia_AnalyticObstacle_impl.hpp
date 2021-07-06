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

        class Fun : public Configurable {
        public:
            utopia::SymbolicFunction gap;
            utopia::SymbolicFunction normal_x;
            utopia::SymbolicFunction normal_y;
            utopia::SymbolicFunction normal_z;

            Fun() : gap("0.1"), normal_x("0"), normal_y("1"), normal_z("0") {}

            void read(Input &in) override {
                std::string expr_gap = "y-0.1";
                std::string expr_normal_x = "0";
                std::string expr_normal_y = "-1";
                std::string expr_normal_z = "0";

                in.get("gap", expr_gap);
                in.get("normal_x", expr_normal_x);
                in.get("normal_y", expr_normal_y);
                in.get("normal_z", expr_normal_z);

                gap = utopia::SymbolicFunction(expr_gap);
                normal_x = utopia::SymbolicFunction(expr_normal_x);
                normal_y = utopia::SymbolicFunction(expr_normal_y);
                normal_z = utopia::SymbolicFunction(expr_normal_z);
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

        {
            Write<Vector> wg(impl_->gap->data()), wn(impl_->normals->data());

            space.node_eval([this, dim](const SizeType idx, const Scalar *point) {
                Scalar coords[4] = {0.0, 0.0, 0.0, 0.0};
                Scalar normal[3] = {0.0, 0.0, 0.0};

                for (int d = 0; d < dim; ++d) {
                    coords[d] = point[d];
                }

                auto g = impl_->fun.gap.eval(coords[0], coords[1], coords[2], coords[3]);

                normal[0] = impl_->fun.normal_x.eval(coords[0], coords[1], coords[2], coords[3]);
                normal[1] = impl_->fun.normal_y.eval(coords[0], coords[1], coords[2], coords[3]);
                normal[2] = impl_->fun.normal_z.eval(coords[0], coords[1], coords[2], coords[3]);

                impl_->gap->data().set(idx * dim, g);

                for (int d = 1; d < dim; ++d) {
                    impl_->gap->data().set(idx * dim + d, 10000000);
                }

                for (int d = 0; d < dim; ++d) {
                    impl_->normals->data().set(idx * dim + d, normal[d]);
                }
            });
        }

        space.apply_constraints(impl_->gap->data());

        impl_->is_contact.zeros(layout(normals->data()));

        int n_var = space.n_var();

        auto r = local_range_device(impl_->is_contact);
        RangeDevice<Vector> rd(r.begin(), r.begin() + r.extent() / n_var);
        auto view = local_view_device(impl_->is_contact);

        parallel_for(
            rd, UTOPIA_LAMBDA(const SizeType i) { view.set(i * n_var, 1.); });

        space.apply_zero_constraints(impl_->is_contact);
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

        space.write("is_contact.e", impl_->is_contact);

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

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
            virtual ~Fun() = default;
            virtual void eval(const Scalar point[3], Scalar &g, Scalar normal[3]) = 0;
            void read(Input &in) override { in.get("side", side); }

            std::string side;
        };

        class Plane : public Fun {
        public:
            using Super = Fun;

            static const std::string name() { return "plane"; }

            Scalar d;
            Scalar normal[3];

            Plane() : d(0) {
                normal[0] = 0;
                normal[1] = 1;
                normal[2] = 0;
            }

            void eval(const Scalar point[3], Scalar &g, Scalar normal[3]) override {
                g = 0.0;

                for (int k = 0; k < 3; ++k) {
                    auto x = point[k] * this->normal[k];
                    g += x;
                }

                g -= d;

                normal[0] = -this->normal[0];
                normal[1] = -this->normal[1];
                normal[2] = -this->normal[2];
            }

            void read(Input &in) override {
                Super::read(in);

                in.get("d", d);

                in.get("normal_x", normal[0]);
                in.get("normal_y", normal[1]);
                in.get("normal_z", normal[2]);
            }
        };

#ifdef UTOPIA_ENABLE_TINY_EXPR
        class Gap : public Fun {
        public:
            using Super = Fun;

            static const std::string name() { return "gap"; }

            utopia::SymbolicFunction gap;

            utopia::SymbolicFunction normal_x;
            utopia::SymbolicFunction normal_y;
            utopia::SymbolicFunction normal_z;

            Gap() : gap("0"), normal_x("0"), normal_y("-1"), normal_z("0") {}

            void eval(const Scalar point[3], Scalar &g, Scalar normal[3]) override {
                g = -gap.eval(point[0], point[1], point[2]);

                normal[0] = normal_x.eval(point[0], point[1], point[2]);
                normal[1] = normal_y.eval(point[0], point[1], point[2]);
                normal[2] = normal_z.eval(point[0], point[1], point[2]);
            }

            void read(Input &in) override {
                Super::read(in);

                std::string expr_gap = "-y";

                std::string expr_normal_x = "0";
                std::string expr_normal_y = "-1";
                std::string expr_normal_z = "0";

                in.get("d", expr_gap);

                in.get("normal_x", expr_normal_x);
                in.get("normal_y", expr_normal_y);
                in.get("normal_z", expr_normal_z);

                gap = utopia::SymbolicFunction(expr_gap);

                normal_x = utopia::SymbolicFunction(expr_normal_x);
                normal_y = utopia::SymbolicFunction(expr_normal_y);
                normal_z = utopia::SymbolicFunction(expr_normal_z);
            }
        };
#endif  // UTOPIA_ENABLE_TINY_EXPR

        std::vector<std::unique_ptr<Fun>> functions;
        bool update_transfer{true};
    };

    template <class FunctionSpace>
    AnalyticObstacle<FunctionSpace>::AnalyticObstacle() : impl_(utopia::make_unique<Impl>()) {}

    template <class FunctionSpace>
    AnalyticObstacle<FunctionSpace>::~AnalyticObstacle() = default;

    template <class FunctionSpace>
    void AnalyticObstacle<FunctionSpace>::read(Input &in) {
        in.get("fun", [this](Input &array_node) {
            array_node.get_all([this](Input &node) {
                std::string type;
                node.get("type", type);
                if (type == Impl::Plane::name()) {
                    auto fun = utopia::make_unique<typename Impl::Plane>();
                    fun->read(node);
                    impl_->functions.push_back(std::move(fun));
#ifdef UTOPIA_ENABLE_TINY_EXPR
                } else if (type == Impl::Gap::name()) {
                    auto fun = utopia::make_unique<typename Impl::Gap>();
                    fun->read(node);
                    impl_->functions.push_back(std::move(fun));
#endif  // UTOPIA_ENABLE_TINY_EXPR
                } else {
                    Utopia::Abort("Unsupported obstacle fun type!");
                }
            });
        });

        in.get("export_tensors", impl_->export_tensors);
    }

    template <class FunctionSpace>
    void AnalyticObstacle<FunctionSpace>::describe(std::ostream &os) const {
        UTOPIA_UNUSED(os);
    }

    template <class FunctionSpace>
    bool AnalyticObstacle<FunctionSpace>::init(const std::shared_ptr<FunctionSpace> &domain) {
        impl_->domain = domain;
        return true;
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

        assert(!gap->empty());
        assert(!normals->empty());

        impl_->gap = gap;
        impl_->normals = normals;
        impl_->is_contact.zeros(layout(normals->data()));

        impl_->gap->data().set(impl_->infinity);

        int n_var = space.n_var();

        {
            auto g_view = view_device(impl_->gap->data());
            auto n_view = view_device(impl_->normals->data());
            auto c_view = view_device(impl_->is_contact);

            Range r = range(impl_->gap->data());
            SizeType r_begin = r.begin() / n_var;
            SizeType r_end = r.end() / n_var;

            for (auto &f : impl_->functions) {
                auto fun = [this, dim, g_view, n_view, c_view, n_var, r_begin, r_end, &f](const SizeType idx,
                                                                                          const Scalar *point) {
                    // if (idx < r_begin || idx >= r_end) return;

                    Scalar point3[3] = {0.0, 0.0, 0.0};
                    Scalar normal3[3] = {0.0, 0.0, 0.0};
                    Scalar g = 0;

                    for (int d = 0; d < dim; ++d) {
                        point3[d] = point[d];
                    }

                    f->eval(point3, g, normal3);
                    g_view.set(idx * n_var, g);
                    c_view.set(idx * n_var, 1);

                    for (int d = 0; d < dim; ++d) {
                        n_view.set(idx * n_var + d, normal3[d]);
                    }
                };

                if (f->side.empty()) {
                    space.node_eval(fun);
                } else {
                    space.node_eval(f->side, fun);
                }
            }
        }

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

        // space.write("gap.e", e_mul(impl_->is_contact, impl_->gap->data()));
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
        out = transpose(impl_->orthogonal_trafo) * in * impl_->orthogonal_trafo;
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

    template <class FunctionSpace>
    std::shared_ptr<typename Traits<FunctionSpace>::Matrix>
    AnalyticObstacle<FunctionSpace>::orthogonal_transformation() {
        return make_ref(impl_->orthogonal_trafo);
    }

}  // namespace utopia

#endif  // UTOPIA_ANALITYC_OBSTACLE_IMPL_HPP

#ifndef UTOPIA_POROUS_FLOW_FE_HPP
#define UTOPIA_POROUS_FLOW_FE_HPP


#include "utopia_Function.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_ProjectionView.hpp"

#include "utopia_QuadratureView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_NodalInterpolateView.hpp"
#include "utopia_Algorithms.hpp"
#include "utopia_DeviceReduce.hpp"
#include "utopia_DeviceTensorReduce.hpp"
#include "utopia_Coefficient.hpp"
#include "utopia_Edge2.hpp"
#include "utopia_ElementWisePseudoInverse.hpp"

namespace utopia {

    template<class FunctionSpace>
    class PorousFlowFE final : public Function<
                            typename FunctionSpace::Matrix,
                            typename FunctionSpace::Vector>,
                            // public Operator<typename FunctionSpace::Vector>,
                            public Configurable {
    public:
        using Comm       = typename FunctionSpace::Comm;
        using Matrix     = typename FunctionSpace::Matrix;
        using Vector     = typename FunctionSpace::Vector;
        using Scalar     = typename Traits<Vector>::Scalar;
        using SizeType   = typename Traits<Vector>::SizeType;
        using Elem       = typename FunctionSpace::Elem;
        using Quadrature = utopia::Quadrature<Elem, 2>;
        using Laplacian  = utopia::Laplacian<FunctionSpace, Quadrature>;
        using ScaledMassMatrix = utopia::ScaledMassMatrix<FunctionSpace, Quadrature>;
        using Coefficient = utopia::Coefficient<FunctionSpace>;

        using Device     = typename FunctionSpace::Device;

        static const int Dim    = Elem::Dim;
        static const int NNodes = Elem::NNodes;
        using ElementMatrix     = utopia::StaticMatrix<Scalar, NNodes, NNodes>;
        using ElementVector     = utopia::StaticVector<Scalar, NNodes>;
        using Point             = typename FunctionSpace::Point;

        using LKernel           = utopia::LaplacianKernel<Scalar>;

        static const int NQPoints = Quadrature::NPoints;


        void read(Input &in) override {
            bool demo = false;

            in.get("demo", demo);
            in.get("rescale_with_spacing", rescale_with_spacing_);

            if(demo) {
                init_demo_fracture_network();

		bool export_permeability = false;

		in.get("export_permeability", export_permeability);

		if(export_permeability) {
		  rename("permeability", *permeability_field_);
		  space_->write("P.vts", *permeability_field_);
		}
            }
        }

        // inline Size size() const
        // {
        //     const SizeType n_dofs = space_->n_dofs();
        //     return {n_dofs, n_dofs};
        // }

        // inline Size local_size() const
        // {
        //     const SizeType n_dofs = space_->n_local_dofs();
        //     return {n_dofs, n_dofs};
        // }

        // inline Comm &comm() override { return space_->comm(); }
        // inline const Comm &comm() const override { return space_->comm(); }

        // bool apply(const Vector &x, Vector &y) const override {
        //     return false;
        // }

        inline bool value(const Vector &, Scalar &) const override
        {
           return false;
        }

        inline bool update(const Vector &x) override {
            return true;
        }

        inline bool gradient(const Vector &x, Vector &g) const override
        {
            return false;
        }

        inline bool hessian(const Vector &x, Matrix &H) const override
        {
            Chrono c;
            c.start();

            if(empty(H)) {
                space_->create_matrix(H);
            } else {
                H *= 0.0;
            }

            PhysicalGradient<FunctionSpace, Quadrature> grad_temp(*space_, quadrature_);
            Differential<FunctionSpace, Quadrature> differential_temp(*space_, quadrature_);
            auto p_val = permeability_field_fun_->value(quadrature_);

           {
               auto space_view = space_->view_device();

               auto dx_view    = differential_temp.view_device();
               auto grad_view  = grad_temp.view_device();

               auto H_view     = space_->assembly_view_device(H);
               auto permeability_view = p_val.view_device();

               Device::parallel_for(
                   space_->local_element_range(),
                   UTOPIA_LAMBDA(const SizeType &i)
               {
                   Elem e;
                   StaticVector<Scalar, NQPoints> permeability;
                   ElementMatrix el_mat;
                   space_view.elem(i, e);
                   permeability_view.get(e, permeability);
                   el_mat.set(0.0);

                   auto grad = grad_view.make(e);
                   auto dx   = dx_view.make(e);

                   const auto n_qp  = grad.n_points();
                   const auto n_fun = grad.n_functions();

                   for(SizeType k = 0; k < n_qp; ++k) {
                       auto ck = permeability(k);

                       for(SizeType j = 0; j < n_fun; ++j) {
                           const auto g_test  = grad(j, k);
                           el_mat(j, j) += LKernel::apply(ck, g_test, g_test, dx(k));

                           for(SizeType l = j + 1; l < n_fun; ++l) {
                               const auto g_trial = grad(l, k);
                               const Scalar v = LKernel::apply(ck, g_trial, g_test, dx(k));

                               el_mat(j, l) += v;
                               el_mat(l, j) += v;
                           }
                       }
                   }

                   space_view.add_matrix(e, el_mat, H_view);
               });
           }

            space_->apply_constraints(H);

            c.stop();
            if(x.comm().rank() == 0) { std::cout << "PoissonFE::hessian(...): " << c << std::endl; }
            return true;
        }

        Vector &permeability_field()
        {
            return *permeability_field_;
        }

        inline bool has_preconditioner() const override
        {
            return false;
        }

        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const
        {
            space_->create_matrix(H);
            return true;
        }

        PorousFlowFE(FunctionSpace &space)
        : space_(utopia::make_ref(space)),
          quadrature_(),
          rescale_with_spacing_(false)
        {
            init();
        }

        PorousFlowFE(const std::shared_ptr<FunctionSpace> &space)
        : space_(space),
          quadrature_(),
          rescale_with_spacing_(false)
        {
            init();
        }

        class FractureNetwork {
        public:
            class LineFracture {
            public:
                Point p1, p2;
                Scalar aperture;
                Scalar permeability;
            };

            void init_demo()
            {
                line_fractures.resize(1);
                auto &p1 = line_fractures[0].p1;
                auto &p2 = line_fractures[0].p2;

                p1[0] = 0.00001;
                p1[1] = 0.2;

                p2[0] = 0.99999;
                p2[1] = 0.8;
            }

            std::vector<LineFracture> line_fractures;
        };

        void init_demo_fracture_network()
        {
            using Point1 = utopia::StaticVector<Scalar, 1>;

            Scalar backround_perm = 1e-4;
            Scalar frac_aperture  = 1e-4;
            Scalar frac_perm      = 1e4;
            Scalar spacing        = space_->mesh().min_spacing();

            ArrayView<Point1, 12> q_points;
            ArrayView<Scalar, 12> q_weights;

            utopia::Quadrature<Scalar, 6, 1>::get(q_points, q_weights);

            Point p1, p2;

            p1[0] = 0.00001;
            p1[1] = 0.2;

            p2[0] = 0.99999;
            p2[1] = 0.8;

            Edge2<Scalar, 2> fracture(p1, p2);

            auto &mesh = space_->mesh();

            mass_vector_ = std::make_shared<Vector>();
            space_->create_vector(*mass_vector_);
            permeability_field_->set(0.0);

            ShapeFunction<FunctionSpace, Quadrature> fun_temp(*space_, quadrature_);
            Differential<FunctionSpace, Quadrature> differential_temp(*space_, quadrature_);
            // PhysicalPoint<FunctionSpace, Quadrature> points_temp(*space_, quadrature_);

            {
                auto space_view = space_->view_device();

                auto dx_view    = differential_temp.view_device();
                auto fun_view   = fun_temp.view_device();

                auto p_view     = space_->assembly_view_device(*permeability_field_);
                auto m_view     = space_->assembly_view_device(*mass_vector_);

                // auto points_view = points_temp.view_device();

                Device::parallel_for(
                    space_->local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                {
                    Elem e;
                    Point p, v, p_quad, isect_1, isect_2;
                    StaticVector<Scalar, 1> p_fracture;
                    StaticVector<Scalar, NQPoints> permeability;
                    ElementVector p_el_vec, m_el_vec;

                    p_el_vec.set(0.0);
                    m_el_vec.set(0.0);

                    space_view.elem(i, e);

                    auto fun    = fun_view.make(e);
                    auto dx     = dx_view.make(e);

                    const auto n_qp  = fun.n_points();
                    const auto n_fun = fun.n_functions();

                    if(e.univar_elem().intersect_line(fracture.node(0), fracture.node(1), isect_1, isect_2)) {
                        UTOPIA_DEVICE_ASSERT(fracture.contains(isect_1));
                        UTOPIA_DEVICE_ASSERT(fracture.contains(isect_2));

                        v = isect_2 - isect_1;

                        for(SizeType k = 0; k < n_qp; ++k) {
                            Scalar w = q_weights[k] * fracture.measure();
                            p = q_points[k](0) * v;
                            fracture.inverse_transform(p, p_fracture);
                            e.inverse_transform(p, p_quad);

                            for(SizeType j = 0; j < n_fun; ++j) {
                                for(SizeType l = 0; l < fracture.n_functions(); ++l) {
                                    const Scalar mm = e.fun(j, p_quad) * fracture.fun(l, p_fracture);
                                    p_el_vec(j) += frac_perm * frac_aperture * mm * w;
                                    m_el_vec(j) += mm * w;
                                }
                            }
                        }

                        space_view.add_vector(e, p_el_vec, p_view);
                        space_view.add_vector(e, m_el_vec, m_view);
                    }
                });
            }

            e_pseudo_inv(*mass_vector_, *mass_vector_, 1e-12);
            (*permeability_field_) = e_mul((*permeability_field_) , (*mass_vector_));
            permeability_field_->shift(backround_perm);
        }

        void init_demo_fracture_network_gaussian()
        {
            Scalar backround_perm = 1e-4;
            Scalar frac_aperture  = 1e-4;
            Scalar frac_perm      = 1e4;
            Scalar spacing        = space_->mesh().min_spacing();

            if(rescale_with_spacing_ && (spacing > frac_aperture)) {
                std::cout << "RESCALING" << std::endl;
                Scalar spacing2 = std::sqrt(frac_aperture/spacing);
                frac_aperture /= spacing2;
                frac_perm *= spacing2;
            }

            std::cout << spacing << std::endl;

            Point p1, p2, u, n;

            p1[0] = 0.00001;
            p1[1] = 0.2;

            p2[0] = 0.99999;
            p2[1] = 0.8;

            // p1[0] = 0.2;
            // p1[1] = 0.2;

            // p2[0] = 0.8;
            // p2[1] = 0.8;

            u = p2 - p1;
            const Scalar len_u = norm2(u);
            u /= len_u;

            n[0] = -u[1];
            n[1] =  u[0];


            // Scalar frac_aperture = 1e-2;


            auto &mesh = space_->mesh();

            // mesh.find_cell(p1);
            // mesh.find_cell(p2);

            mass_vector_ = std::make_shared<Vector>();
            space_->create_vector(*mass_vector_);
            permeability_field_->set(0.0);

            ShapeFunction<FunctionSpace, Quadrature> fun_temp(*space_, quadrature_);
            Differential<FunctionSpace, Quadrature> differential_temp(*space_, quadrature_);
            PhysicalPoint<FunctionSpace, Quadrature> points_temp(*space_, quadrature_);

            {
                auto space_view = space_->view_device();

                auto dx_view    = differential_temp.view_device();
                auto fun_view   = fun_temp.view_device();

                auto p_view     = space_->assembly_view_device(*permeability_field_);
                auto m_view     = space_->assembly_view_device(*mass_vector_);

                auto points_view = points_temp.view_device();

                Device::parallel_for(
                    space_->local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                {
                    Elem e;
                    Point p, v;
                    StaticVector<Scalar, NQPoints> permeability;
                    ElementVector p_el_vec, m_el_vec;

                    p_el_vec.set(0.0);
                    m_el_vec.set(0.0);

                    space_view.elem(i, e);

                    auto fun    = fun_view.make(e);
                    auto dx     = dx_view.make(e);
                    auto points = points_view.make(e);

                    const auto n_qp  = fun.n_points();
                    const auto n_fun = fun.n_functions();

                    for(SizeType k = 0; k < n_qp; ++k) {
                        points.get(k, p);
                        v = p - p1;

                        Scalar v_t = dot(u, v);
                        if(v_t < 0) {
                            v_t = -v_t;
                        } else if(v_t > 0) {
                            if(v_t <= len_u) {
                                v_t = 0.0;
                            } else {
                                v_t -= len_u;
                            }
                        }

                        const Scalar v_n   = dot(n, v);
                        const Scalar dist  = device::sqrt(v_t * v_t + v_n * v_n);
                        const Scalar weight = device::exp(-0.5 * (dist*dist)/(frac_aperture * frac_aperture));
                        const Scalar perm_k = weight * frac_perm + (1-weight) * backround_perm;

                        for(SizeType j = 0; j < n_fun; ++j) {
                            p_el_vec(j) += perm_k * fun(j, k) * dx(k);
                            m_el_vec(j) += fun(j, k) * dx(k);
                        }
                    }

                    space_view.add_vector(e, p_el_vec, p_view);
                    space_view.add_vector(e, m_el_vec, m_view);
                });
            }

            const Scalar avg_perm = sum(*permeability_field_);
            std::cout << "avg_perm: " << avg_perm << std::endl;
            (*permeability_field_) = e_mul((*permeability_field_) , 1./(*mass_vector_));
        }

    private:
        std::shared_ptr<FunctionSpace> space_;
        Quadrature quadrature_;

        std::shared_ptr<Vector> permeability_field_;
        std::shared_ptr<Vector> mass_vector_;
        std::unique_ptr<FEFunction<FunctionSpace> > permeability_field_fun_;
        bool rescale_with_spacing_;

        void init()
        {
            permeability_field_ = std::make_shared<Vector>();
            space_->create_vector(*permeability_field_);
            permeability_field_->set(1.0);
            permeability_field_fun_ = utopia::make_unique<FEFunction<FunctionSpace>>(space_, permeability_field_);
        }
    };

}



#endif //UTOPIA_POROUS_FLOW_FE_HPP

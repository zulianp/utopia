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
            auto p_val = porosity_field_fun_->value(quadrature_);

           {
               auto space_view = space_->view_device();

               auto dx_view    = differential_temp.view_device();
               auto grad_view  = grad_temp.view_device();

               auto H_view     = space_->assembly_view_device(H);
               auto porosity_view = p_val.view_device();

               Device::parallel_for(
                   space_->local_element_range(),
                   UTOPIA_LAMBDA(const SizeType &i)
               {
                   Elem e;
                   StaticVector<Scalar, NQPoints> porosity;
                   ElementMatrix el_mat;
                   space_view.elem(i, e);
                   porosity_view.get(e, porosity);
                   el_mat.set(0.0);

                   auto grad = grad_view.make(e);
                   auto dx   = dx_view.make(e);

                   const auto n_qp  = grad.n_points();
                   const auto n_fun = grad.n_functions();

                   for(SizeType k = 0; k < n_qp; ++k) {
                       auto ck = porosity(k);

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

        Vector &porosity_field()
        {
            return *porosity_field_;
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
          quadrature_()
        {
            init();
        }

        PorousFlowFE(const std::shared_ptr<FunctionSpace> &space)
        : space_(space),
          quadrature_()
        {
            init();
        }

    private:
        std::shared_ptr<FunctionSpace> space_;
        Quadrature quadrature_;

        std::shared_ptr<Vector> porosity_field_;
        std::shared_ptr<Vector> mass_vector_;
        std::unique_ptr<FEFunction<FunctionSpace> > porosity_field_fun_;

        void init_demo_fracture_network()
        {
            Scalar backround_perm = 1e-4;

            Point p1, p2;

            p1[0] = 0.0;
            p1[1] = 0.2;

            p2[0] = 1.0;
            p2[1] = 0.9;

            Scalar frac_aperture = 1e-4;
            Scalar frac_perm     = 1e4;

            auto &mesh = space_->mesh();

            // mesh.find_cell(p1);
            // mesh.find_cell(p2);

            mass_vector_ = std::make_shared<Vector>();
            space_->create_vector(*mass_vector_);


            ShapeFunction<FunctionSpace, Quadrature> fun_temp(*space_, quadrature_);
            Differential<FunctionSpace, Quadrature> differential_temp(*space_, quadrature_);

            {
                auto space_view = space_->view_device();

                auto dx_view    = differential_temp.view_device();
                auto fun_view   = fun_temp.view_device();

                auto p_view     = space_->assembly_view_device(*porosity_field_);
                auto m_view     = space_->assembly_view_device(*mass_vector_);

                Device::parallel_for(
                    space_->local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                {
                    Elem e;
                    StaticVector<Scalar, NQPoints> porosity;
                    ElementVector p_el_vec, m_el_vec;

                    p_el_vec.set(0.0);
                    m_el_vec.set(0.0);

                    space_view.elem(i, e);

                    auto fun  = fun_view.make(e);
                    auto dx   = dx_view.make(e);

                    const auto n_qp  = fun.n_points();
                    const auto n_fun = fun.n_functions();

                    // for(SizeType k = 0; k < n_qp; ++k) {
                    //     auto ck = porosity(k);

                    //     for(SizeType j = 0; j < n_fun; ++j) {
                    //         const auto g_test  = grad(j, k);
                    //         el_mat(j, j) += LKernel::apply(ck, g_test, g_test, dx(k));

                    //         for(SizeType l = j + 1; l < n_fun; ++l) {
                    //             const auto g_trial = grad(l, k);
                    //             const Scalar v = LKernel::apply(ck, g_trial, g_test, dx(k));

                    //             el_mat(j, l) += v;
                    //             el_mat(l, j) += v;
                    //         }
                    //     }
                    // }

                    space_view.add_vector(e, p_el_vec, p_view);
                    space_view.add_vector(e, m_el_vec, m_view);
                });
            }
        }

        void init()
        {
            porosity_field_ = std::make_shared<Vector>();
            space_->create_vector(*porosity_field_);
            porosity_field_->set(1.0);
            porosity_field_fun_ = utopia::make_unique<FEFunction<FunctionSpace>>(space_, porosity_field_);
        }
    };

}



#endif //UTOPIA_POROUS_FLOW_FE_HPP

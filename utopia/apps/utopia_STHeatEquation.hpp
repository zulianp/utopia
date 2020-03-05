#ifndef UTOPIA_ST_HEAT_EQUATION_HPP
#define UTOPIA_ST_HEAT_EQUATION_HPP

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
#include "utopia_SpaceTimeDeriv.hpp"

namespace utopia {

    template<class FunctionSpace>
    class STHeatEquation final : public Function<
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
        static const int NFunctions = Elem::NFunctions;
        static const int NQuadPoints = Quadrature::NPoints;


        using ElementMatrix     = utopia::StaticMatrix<Scalar, NNodes, NNodes>;
        using ElementVector     = utopia::StaticVector<Scalar, NNodes>;
        using Point             = typename FunctionSpace::Point;

        using LKernel           = utopia::LaplacianKernel<Scalar>;


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

            SpaceTimeDeriv<FunctionSpace, Quadrature> st_deriv(*space_, quadrature_);
            // ShapeFunction<FunctionSpace, Quadrature> shape_fun(*space_, quadrature_);
            Differential<FunctionSpace, Quadrature> differential(*space_, quadrature_);

            {
                auto space_view     = space_->view_device();
                auto H_view         = space_->assembly_view_device(H);
                auto st_deriv_view  = st_deriv.view_device();
                // auto shape_fun_view = shape_fun.view_device();
                auto dx_view        = differential.view_device();

                Device::parallel_for(
                    space_->local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                {
                    ElementMatrix el_mat; el_mat.set(0.0);

                    Elem e;
                    space_view.elem(i, e);

                    auto &&dx    = dx_view.make(e);
                    auto &&deriv = st_deriv_view.make(e);
                    // auto &&shape = shape_fun_view.make(e);

                    for(SizeType qp = 0; qp < NQuadPoints; ++qp) {
                        for(SizeType j = 0; j < NFunctions; ++j) {
                            for(SizeType l = 0; l < NFunctions; ++l) {

                                el_mat(j, l) += (
                                    deriv.partial_t(l, qp) * deriv.partial_t(j, qp)
                                    + inner(deriv.grad_x(l, qp), deriv.grad_x_partial_t(j, qp))
                                ) * dx(qp);
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

        inline bool has_preconditioner() const override
        {
            return false;
        }

        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const
        {
            space_->create_matrix(H);
            return true;
        }




        STHeatEquation(FunctionSpace &space)
        : space_(utopia::make_ref(space)),
          quadrature_()
        {}

        STHeatEquation(const std::shared_ptr<FunctionSpace> &space)
        : space_(space),
          quadrature_()
        {}

    private:
        std::shared_ptr<FunctionSpace> space_;
        Quadrature quadrature_;
    };

}


#endif
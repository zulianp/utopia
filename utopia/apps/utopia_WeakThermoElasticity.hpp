#ifndef UTOPIA_WEAK_THERMO_ELASTICITY_HPP
#define UTOPIA_WEAK_THERMO_ELASTICITY_HPP

#include "utopia_LinearElasticityView.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_Views.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_DiffController.hpp"
#include "utopia_StrainView.hpp"

namespace utopia {

    template<class FunctionSpace, int Dim = FunctionSpace::Dim>
    class WeakThermoElasticity final : public Function<
                            typename FunctionSpace::Matrix,
                            typename FunctionSpace::Vector>,
                            public Configurable {
    public:
        using Scalar    = typename FunctionSpace::Scalar;
        using SizeType  = typename FunctionSpace::SizeType;
        using Vector    = typename FunctionSpace::Vector;
        using Matrix    = typename FunctionSpace::Matrix;
        using Device    = typename FunctionSpace::Device;

        using USpace    = typename FunctionSpace::template Subspace<Dim>;
        using CSpace    = typename FunctionSpace::template Subspace<1>;

        using UElem     = typename USpace::ViewDevice::Elem;
        using CElem     = typename CSpace::ViewDevice::Elem;
        using MixedElem = typename FunctionSpace::ViewDevice::Elem;

        //FIXME
        using Quadrature = utopia::Quadrature<typename FunctionSpace::Shape, 2>;

        using CLaplacian  = utopia::Laplacian<CSpace, Quadrature>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        using LEKernel    = utopia::LinearElasticityKernel<Scalar>;


        WeakThermoElasticity(FunctionSpace &space)
        : space_(space)
        {}

        void read(Input &in) override
        {

        }

        inline bool update(const Vector &x) override {
            // x_coeff_.update(x);
            return true;
        }

        bool value(
            const Vector &x_const,
            Scalar &val
        ) const override
        {
            assert(false);
            return false;
        }

        bool gradient(const Vector &x_const, Vector &g) const override
        {
            return false;
        }

        bool hessian(const Vector &x_const, Matrix &H) const override
        {
            if(empty(H)) {
                space_.create_matrix(H);
            } else {
                H *= 0.0;
            }

            USpace U;
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

            auto &x = const_cast<Vector &>(x_const);

            FEFunction<CSpace> c_fun(C, x);
            FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto c_val  = c_fun.value(q);
            auto c_grad = c_fun.gradient(q);
            auto u_val  = u_fun.value(q);
            auto differential = C.differential(q);

            auto v_grad_shape = U.shape_grad(q);
            auto c_shape      = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            Strain<USpace, Quadrature> strain(U, q);

            CLaplacian laplacian(C, q);

            const Scalar E = 50e3;
            const Scalar nu = 0.2;
            const Scalar mu = E/2/(1+nu);
            const Scalar lmbda = E*nu/(1+nu)/(1-2*nu);
            const Scalar alpha = 1e-5;

            {
                auto U_view      = U.view_device();
                auto C_view      = C.view_device();
                auto space_view  = space_.view_device();

                auto c_view      = c_val.view_device();
                auto c_grad_view = c_grad.view_device();
                auto u_view      = u_val.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto v_grad_shape_view = v_grad_shape.view_device();
                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto H_view = space_.assembly_view_device(H);
                auto laplacian_view = laplacian.view_device();

                Device::parallel_for(
                    space_.local_element_range(),
                    UTOPIA_LAMBDA(const SizeType &i)
                    {
                        StaticMatrix<Scalar, C_NDofs, C_NDofs> lap_mat;
                        StaticMatrix<Scalar, U_NDofs, U_NDofs> elast_mat;
                        StaticMatrix<Scalar, U_NDofs, C_NDofs> coupling_mat;
                        StaticMatrix<Scalar, U_NDofs + C_NDofs, U_NDofs + C_NDofs> el_mat;

                        lap_mat.set(0.0);
                        elast_mat.set(0.0);
                        coupling_mat.set(0.0);
                        el_mat.set(0.0);

                        ////////////////////////////////////////////////////////////////////

                        MixedElem e;
                        space_view.elem(i, e);

                        CElem c_e;
                        C_view.elem(i, c_e);

                        auto dx = differential_view.make(c_e);
                        auto shape_c = c_shape_view.make(c_e);

                        UElem u_e;
                        U_view.elem(i, u_e);

                        //(0, 0)
                        laplacian_view.assemble(i, c_e, lap_mat);

                        auto &&strain = strain_view.make(u_e);

                        for(SizeType qp = 0; qp < NQuadPoints; ++qp) {
                            const auto dx_qp = dx(qp);

                            for(SizeType i = 0; i < U_NDofs; ++i) {
                                auto &&strain_test = strain(i, qp);

                                elast_mat(i, i) +=
                                    LEKernel::strain_apply(
                                        mu,
                                        lmbda,
                                        strain_test,
                                        strain_test,
                                        dx_qp
                                );


                                //a(u, u)
                                for(SizeType j = i + 1; j < U_NDofs; ++j) {
                                    auto &&strain_trial = strain(j, qp);

                                    auto v = LEKernel::strain_apply(
                                        mu,
                                        lmbda,
                                        strain_trial,
                                        strain_test,
                                        dx_qp
                                    );

                                    elast_mat(i, j) += v;
                                    elast_mat(j, i) += v;
                                }

                                //a(u, c)
                                for(SizeType j = 0; j < C_NDofs; ++j) {
                                    auto v = -alpha*(3*lmbda+2*mu) * shape_c(j, qp);
                                    coupling_mat(i, j) += inner(strain_test, v * device::identity<Scalar>()) * dx_qp;
                                }
                            }
                        }


                        //to temporary mat buffer
                        el_mat.set_matrix(0, 0, lap_mat);
                        el_mat.set_matrix(C_NDofs, C_NDofs, elast_mat);
                        //In the weak coupling this is not here
                        // el_mat.set_matrix(0, C_NDofs, transpose(coupling_mat));
                        el_mat.set_matrix(C_NDofs, 0, coupling_mat);

                        space_view.add_matrix(e, el_mat, H_view);
                    }
                );
            }

            space_.apply_constraints(H);
            return false;
        }

    private:
        FunctionSpace &space_;
    };

}

#endif
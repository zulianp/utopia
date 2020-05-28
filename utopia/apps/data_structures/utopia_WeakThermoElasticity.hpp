#ifndef UTOPIA_WEAK_THERMO_ELASTICITY_HPP
#define UTOPIA_WEAK_THERMO_ELASTICITY_HPP

#include "utopia_DiffController.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_PrincipalShapeStressView.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_Views.hpp"

namespace utopia {

    // https://comet-fenics.readthedocs.io/en/latest/demo/thermoelasticity/thermoelasticity.html

    template <class FunctionSpace, int Dim = FunctionSpace::Dim>
    class WeakThermoElasticity final : public Function<typename FunctionSpace::Matrix, typename FunctionSpace::Vector>,
                                       public Operator<typename FunctionSpace::Vector> {
    public:
        using Comm = typename FunctionSpace::Comm;
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Device = typename FunctionSpace::Device;

        using USpace = typename FunctionSpace::template Subspace<Dim>;
        using CSpace = typename FunctionSpace::template Subspace<1>;

        using UElem = typename USpace::ViewDevice::Elem;
        using CElem = typename CSpace::ViewDevice::Elem;
        using MixedElem = typename FunctionSpace::ViewDevice::Elem;
        using Coefficient = utopia::Coefficient<FunctionSpace>;

        // FIXME
        using Quadrature = utopia::Quadrature<typename FunctionSpace::Shape, 2>;

        using CLaplacian = utopia::Laplacian<CSpace, Quadrature>;

        static const int C_NDofs = CSpace::NDofs;
        static const int U_NDofs = USpace::NDofs;

        static const int NQuadPoints = Quadrature::NPoints;

        using LEKernel = utopia::LinearElasticityKernel<Scalar>;

        WeakThermoElasticity(FunctionSpace &space) : space_(space), x_coeff_(utopia::make_unique<Coefficient>(space_)) {
            const Scalar E = 50e3;
            const Scalar nu = 0.2;
            mu_ = E / 2 / (1 + nu);
            lambda_ = E * nu / (1 + nu) / (1 - 2 * nu);
            alpha_ = 1e-5;

            H_ptr_ = utopia::make_unique<Matrix>();
        }

        void read(Input &) override {}

        inline bool update(const Vector & /*x*/) override {
            // x_coeff_.update(x);
            return true;
        }

        bool value(const Vector & /*x_const*/, Scalar & /*val*/
                   ) const override {
            assert(false);
            return false;
        }

        bool gradient(const Vector & /*x_const*/, Vector & /*g*/) const override { return false; }

        inline Size size() const override {
            const SizeType n_dofs = space_.n_dofs();
            return {n_dofs, n_dofs};
        }

        inline Size local_size() const override {
            const SizeType n_dofs = space_.n_local_dofs();
            return {n_dofs, n_dofs};
        }

        inline Comm &comm() override { return space_.comm(); }
        inline const Comm &comm() const override { return space_.comm(); }

        bool apply(const Vector &x, Vector &y) const override {
            // const auto &comm = space_.comm();

            if (y.empty()) {
                space_.create_vector(y);
            } else {
                y.set(0.0);
            }

            x_coeff_->update(x);

            USpace U;
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

            Quadrature q;

            auto differential = C.differential(q);

            auto c_shape = C.shape(q);
            auto c_grad_shape = C.shape_grad(q);

            Strain<USpace, Quadrature> strain(U, q);

            CLaplacian laplacian(C, q);

            const Scalar mu = mu_;
            const Scalar lmbda = lambda_;
            const Scalar alpha = alpha_;

            {
                auto x_view = x_coeff_->view_device();
                auto y_view = space_.assembly_view_device(y);

                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto space_view = space_.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto c_shape_view = c_shape.view_device();
                auto c_grad_shape_view = c_grad_shape.view_device();

                auto g_view = space_.assembly_view_device(y);
                // auto laplacian_view = laplacian.view_device();

                Device::parallel_for(space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                    StaticVector<Scalar, C_NDofs + U_NDofs> coeff, el_vec;
                    el_vec.set(0.0);

                    ////////////////////////////////////////////////////////////////////

                    MixedElem e;
                    space_view.elem(i, e);

                    x_view.get(e, coeff);

                    CElem c_e;
                    C_view.elem(i, c_e);

                    auto dx = differential_view.make(c_e);
                    auto shape_c = c_shape_view.make(c_e);
                    auto grad_c = c_grad_shape_view.make(c_e);

                    UElem u_e;
                    U_view.elem(i, u_e);

                    auto &&strain = strain_view.make(u_e);

                    for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                        const auto dx_qp = dx(qp);

                        // laplacian
                        for (SizeType i = 0; i < C_NDofs; ++i) {
                            auto &&g_test = grad_c(i, qp);

                            el_vec(i) += inner(g_test, g_test) * dx_qp * coeff(i);

                            for (SizeType j = i + 1; j < C_NDofs; ++j) {
                                auto &&g_trial = grad_c(j, qp);

                                auto v = inner(g_test, g_trial) * dx_qp;

                                el_vec(i) += v * coeff(j);
                                el_vec(j) += v * coeff(i);
                            }
                        }

                        for (SizeType i = 0; i < U_NDofs; ++i) {
                            auto &&strain_test = strain(i, qp);

                            el_vec(C_NDofs + i) +=
                                LEKernel::strain_apply(mu, lmbda, strain_test, strain_test, dx_qp) * coeff(C_NDofs + i);

                            // a(u, u)
                            for (SizeType j = i + 1; j < U_NDofs; ++j) {
                                auto &&strain_trial = strain(j, qp);

                                auto v = LEKernel::strain_apply(mu, lmbda, strain_trial, strain_test, dx_qp);

                                el_vec(C_NDofs + i) += v * coeff(C_NDofs + j);
                                el_vec(C_NDofs + j) += v * coeff(C_NDofs + i);
                            }

                            // a(u, c)
                            for (SizeType j = 0; j < C_NDofs; ++j) {
                                auto v = -alpha * (3 * lmbda + 2 * mu) * shape_c(j, qp);
                                el_vec(C_NDofs + i) +=
                                    inner(strain_test, v * device::identity<Scalar>()) * dx_qp * coeff(j);
                            }
                        }
                    }

                    space_view.add_vector(e, el_vec, y_view);
                });
            }

            space_.copy_at_constrained_dofs(x, y);

            assert(check_with_matrix(x, y));
            return true;
        }

        bool check_with_matrix(const Vector &x, const Vector &y) const {
            if (H_ptr_->empty()) {
                space_.create_matrix(*H_ptr_);
                hessian(x, *H_ptr_);
            }

            Vector mat_y = (*H_ptr_) * x;

            Scalar norm_diff = norm_infty(y - mat_y);

            bool ok = norm_diff < 1e-8;

            if (!ok) {
                Scalar norm_y = norm_infty(mat_y);

                std::cout << "norm_diff: " << norm_diff << std::endl;
                std::cout << "norm_y: " << norm_y << std::endl;
                // Vector diff = y - mat_y;
                // disp(diff);
            }
            // std::cout << "norm_diff: " << norm_diff << std::endl;
            // assert(norm_diff < 1e-8);
            // return ok;

            return true;
        }

        bool hessian(const Vector &, Matrix &H) const override {
            if (empty(H)) {
                space_.create_matrix(H);
            } else {
                H *= 0.0;
            }

            USpace U;
            space_.subspace(1, U);
            CSpace C = space_.subspace(0);

            // auto &x = const_cast<Vector &>(x_const);
            //
            // FEFunction<CSpace> c_fun(C, x);
            // FEFunction<USpace> u_fun(U, x);

            Quadrature q;

            auto differential = C.differential(q);

            auto c_shape = C.shape(q);
            // auto c_grad_shape = C.shape_grad(q);

            Strain<USpace, Quadrature> strain(U, q);

            CLaplacian laplacian(C, q);

            const Scalar mu = mu_;
            const Scalar lmbda = lambda_;
            const Scalar alpha = alpha_;

            {
                auto U_view = U.view_device();
                auto C_view = C.view_device();
                auto space_view = space_.view_device();

                auto strain_view = strain.view_device();
                auto differential_view = differential.view_device();

                auto c_shape_view = c_shape.view_device();
                // auto c_grad_shape_view = c_grad_shape.view_device();

                auto H_view = space_.assembly_view_device(H);
                auto laplacian_view = laplacian.view_device();

                Device::parallel_for(space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
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

                    for (SizeType qp = 0; qp < NQuadPoints; ++qp) {
                        const auto dx_qp = dx(qp);

                        for (SizeType i = 0; i < U_NDofs; ++i) {
                            auto &&strain_test = strain(i, qp);

                            elast_mat(i, i) += LEKernel::strain_apply(mu, lmbda, strain_test, strain_test, dx_qp);

                            // a(u, u)
                            for (SizeType j = i + 1; j < U_NDofs; ++j) {
                                auto &&strain_trial = strain(j, qp);

                                auto v = LEKernel::strain_apply(mu, lmbda, strain_trial, strain_test, dx_qp);

                                elast_mat(i, j) += v;
                                elast_mat(j, i) += v;
                            }

                            // a(u, c)
                            for (SizeType j = 0; j < C_NDofs; ++j) {
                                auto v = -alpha * (3 * lmbda + 2 * mu) * shape_c(j, qp);
                                coupling_mat(i, j) += inner(strain_test, v * device::identity<Scalar>()) * dx_qp;
                            }
                        }
                    }

                    // to temporary mat buffer
                    el_mat.set_matrix(0, 0, lap_mat);
                    el_mat.set_matrix(C_NDofs, C_NDofs, elast_mat);
                    // In the weak coupling this is not here
                    // el_mat.set_matrix(0, C_NDofs, transpose(coupling_mat));
                    el_mat.set_matrix(C_NDofs, 0, coupling_mat);

                    space_view.add_matrix(e, el_mat, H_view);
                });
            }

            space_.apply_constraints(H);
            return false;
        }

    private:
        FunctionSpace &space_;
        Scalar mu_, lambda_, alpha_;

        ///
        std::unique_ptr<Coefficient> x_coeff_;

        // debugging
        std::unique_ptr<Matrix> H_ptr_;
    };

}  // namespace utopia

#endif
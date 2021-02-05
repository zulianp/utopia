#ifndef UTOPIA_POISSON_FE_HPP
#define UTOPIA_POISSON_FE_HPP

#include "utopia_Function.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_ProjectionView.hpp"

#include "utopia_Algorithms.hpp"
#include "utopia_Coefficient.hpp"
#include "utopia_DeviceReduce.hpp"
#include "utopia_DeviceTensorReduce.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_NodalInterpolateView.hpp"
#include "utopia_QuadratureView.hpp"

#ifdef UTOPIA_WITH_VC
#define USE_SIMD_ASSEMBLY
#endif

#ifdef USE_SIMD_ASSEMBLY
#include "utopia_simd_Assembler.hpp"
#endif  // USE_SIMD_ASSEMBLY

namespace utopia {

#ifndef USE_SIMD_ASSEMBLY
    namespace simd {
        template <typename T>
        inline T integrate(const T v) {
            return v;
        }
    }   // namespace simd
#endif  // USE_SIMD_ASSEMBLY

    template <class FunctionSpace>
    class PoissonFE final : public Function<typename FunctionSpace::Matrix, typename FunctionSpace::Vector>,
                            public Operator<typename FunctionSpace::Vector> {
    public:
        using Comm = typename FunctionSpace::Comm;
        using Matrix = typename FunctionSpace::Matrix;
        using Vector = typename FunctionSpace::Vector;
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Elem = typename FunctionSpace::Elem;

#ifdef USE_SIMD_ASSEMBLY
        using SIMDType = Vc::Vector<Scalar>;
        using Quadrature = simd::Quadrature<SIMDType, Elem::Dim>;
        using GradValue = typename simd::FETraits<Elem, SIMDType>::GradValue;
#else
        using Quadrature = utopia::Quadrature<Elem, 2 * (Elem::Order - 1)>;
        using GradValue = typename Elem::GradValue;
#endif  // USE_SIMD_ASSEMBLY

        using Laplacian = utopia::Laplacian<FunctionSpace, Quadrature>;
        using ScaledMassMatrix = utopia::ScaledMassMatrix<FunctionSpace, Quadrature>;
        using Coefficient = utopia::Coefficient<FunctionSpace>;

        using Device = typename FunctionSpace::Device;

        static const int Dim = Elem::Dim;
        static const int NNodes = Elem::NNodes;
        static const int NFunctions = Elem::NFunctions;
        using ElementMatrix = utopia::StaticMatrix<Scalar, NFunctions, NFunctions>;
        using ElementVector = utopia::StaticVector<Scalar, NFunctions>;
        using Point = typename FunctionSpace::Point;

        using LKernel = utopia::LaplacianKernel<Scalar>;

        void read(Input & /*in*/) override {}

        inline Size size() const override {
            const SizeType n_dofs = space_->n_dofs();
            return {n_dofs, n_dofs};
        }

        inline Size local_size() const override {
            const SizeType n_dofs = space_->n_local_dofs();
            return {n_dofs, n_dofs};
        }

        inline Comm &comm() override { return space_->comm(); }
        inline const Comm &comm() const override { return space_->comm(); }

        bool apply(const Vector &x, Vector &y) const override {
            UTOPIA_TRACE_REGION_BEGIN("PoissonFE::apply(...)");

            if (empty(y)) {
                space_->create_vector(y);
            } else {
                y.set(0.0);
            }

            x_coeff_->update(x);

            PhysicalGradient<FunctionSpace, Quadrature> grad_temp(*space_, quadrature_);
            Differential<FunctionSpace, Quadrature> differential_temp(*space_, quadrature_);

            {
                auto space_view = space_->view_device();
                auto y_view = space_->assembly_view_device(y);
                auto coeff_view = x_coeff_->view_device();

                auto dx_view = differential_temp.view_device();
                auto grad_view = grad_temp.view_device();

                Device::parallel_for(
                    space_->element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        GradValue g_test, g_trial;
                        Elem e;
                        ElementVector coeff, el_vec;
                        space_view.elem(i, e);

                        coeff_view.get(e, coeff);
                        el_vec.set(0.0);

                        auto grad = grad_view.make(e);
                        auto dx = dx_view.make(e);

                        const int n_qp = grad.n_points();
                        const int n_fun = grad.n_functions();

                        for (int k = 0; k < n_qp; ++k) {
                            for (int j = 0; j < n_fun; ++j) {
                                grad.get(j, k, g_test);
                                el_vec(j) += simd::integrate(LKernel::apply(1.0, g_test, g_test, dx(k)) * coeff(j));

                                for (int l = j + 1; l < n_fun; ++l) {
                                    grad.get(l, k, g_trial);
                                    const Scalar v = simd::integrate(LKernel::apply(1.0, g_trial, g_test, dx(k)));

                                    el_vec(j) += v * coeff(l);
                                    el_vec(l) += v * coeff(j);
                                }
                            }
                        }

                        space_view.add_vector(e, el_vec, y_view);
                    });
            }

            space_->copy_at_constrained_dofs(x, y);

            UTOPIA_TRACE_REGION_END("PoissonFE::apply(...)");
            return false;
        }

        inline bool value(const Vector &, Scalar &) const override { return false; }

        inline bool update(const Vector &x) override {
            x_coeff_->update(x);
            return true;
        }

        inline bool gradient(const Vector &x, Vector &g) const override {
            Chrono c;
            c.start();

            if (empty(g)) {
                space_->create_vector(g);
            } else {
                g.set(0.0);
            }

            {
                auto space_view = space_->view_device();

                // auto x_view = space_->assembly_view_device(x);
                auto g_view = space_->assembly_view_device(g);

                auto l_view = laplacian_->view_device();
                auto coeff_view = x_coeff_->view_device();

                Device::parallel_for(
                    space_->element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        Elem e;
                        space_view.elem(i, e);

                        ElementVector coeff;

                        coeff_view.get(e, coeff);

                        ElementMatrix el_mat;
                        el_mat.set(0.0);

                        l_view.assemble(i, e, el_mat);

                        ElementVector el_vec;
                        el_vec = el_mat * coeff;

                        space_view.add_vector(e, el_vec, g_view);
                    });
            }

            if (!empty(rhs_)) {
                g -= rhs_;
            }

            space_->apply_zero_constraints(g);

            c.stop();
            if (x.comm().rank() == 0) {
                utopia::out() << "PoissonFE::gradient(...): " << c << std::endl;
            }
            return true;
        }

        inline bool hessian(const Vector &x, Matrix &H) const override {
            Chrono c;
            c.start();

            if (empty(H)) {
                space_->create_matrix(H);
            } else {
                H *= 0.0;
            }

            {
                auto space_view = space_->view_device();

                auto H_view = space_->assembly_view_device(H);
                auto l_view = laplacian_->view_device();

                Device::parallel_for(
                    space_->element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        Elem e;
                        space_view.elem(i, e);

                        ElementMatrix el_mat;
                        el_mat.set(0.0);
                        l_view.assemble(i, e, el_mat);

                        space_view.add_matrix(e, el_mat, H_view);
                    });
            }

            space_->apply_constraints(H);

            c.stop();
            if (x.comm().rank() == 0) {
                utopia::out() << "PoissonFE::hessian(...): " << c << std::endl;
            }
            return true;
        }

        inline bool has_preconditioner() const override { return false; }

        inline bool initialize_hessian(Matrix &H, Matrix & /*H_pre*/) const override {
            space_->create_matrix(H);
            return true;
        }

        template <class Fun>
        void init_forcing_function(Fun fun) {
            space_->create_vector(rhs_);
            Projection<FunctionSpace, Quadrature, Fun> proj(*space_, quadrature_, fun);

            auto proj_view = proj.view_device();

            {
                auto space_view = space_->view_device();
                auto rhs_view = space_->assembly_view_device(rhs_);

                Device::parallel_for(
                    space_->element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        Elem e;
                        space_view.elem(i, e);

                        ElementVector el_vec;
                        el_vec.set(0.0);

                        proj_view.assemble(i, e, el_vec);
                        space_view.add_vector(e, el_vec, rhs_view);
                    });
            }

            space_->apply_constraints(rhs_);
        }

        PoissonFE(FunctionSpace &space)
            : quadrature_(), space_(utopia::make_ref(space)), laplacian_(nullptr), x_coeff_(nullptr) {
            init();
        }

        PoissonFE(const std::shared_ptr<FunctionSpace> &space)
            : quadrature_(), space_(space), laplacian_(nullptr), x_coeff_(nullptr) {
            init();
        }

    private:
        Quadrature quadrature_;
        std::shared_ptr<FunctionSpace> space_;
        std::unique_ptr<Laplacian> laplacian_;
        std::unique_ptr<Coefficient> x_coeff_;
        Vector rhs_;

        void init() {
            x_coeff_ = utopia::make_unique<Coefficient>(*space_);

#ifdef USE_SIMD_ASSEMBLY
            simd::QuadratureDB<Elem, SIMDType>::get(2 * (Elem::Order - 1), quadrature_);
#else
            laplacian_ = utopia::make_unique<Laplacian>(*space_, quadrature_);

#endif  // USE_SIMD_ASSEMBLY
        }
    };
}  // namespace utopia

#endif  // UTOPIA_POISSON_FE_HPP

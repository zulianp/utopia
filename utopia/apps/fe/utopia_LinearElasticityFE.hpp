#ifndef UTOPIA_LINEAR_ELASTICITY_FE_HPP
#define UTOPIA_LINEAR_ELASTICITY_FE_HPP

#include "utopia_AssemblyView.hpp"
#include "utopia_BratuFE.hpp"
#include "utopia_ConjugateGradient.hpp"
#include "utopia_Core.hpp"
#include "utopia_DeviceView.hpp"
#include "utopia_FEFunction.hpp"
#include "utopia_GradInterpolate.hpp"
#include "utopia_Input.hpp"
#include "utopia_LaplacianView.hpp"
#include "utopia_LinearElasticityView.hpp"
#include "utopia_MPITimeStatistics.hpp"
#include "utopia_MassMatrixView.hpp"
#include "utopia_PhaseField.hpp"
#include "utopia_PoissonFE.hpp"
#include "utopia_PrincipalStrainsView.hpp"
#include "utopia_SampleView.hpp"
#include "utopia_StrainView.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_make_unique.hpp"

// petsc
#include "utopia_petsc.hpp"
#include "utopia_petsc_DM.hpp"
#include "utopia_petsc_DMDA_FunctionSpace.hpp"
#include "utopia_petsc_DirichletBoundaryConditions.hpp"
#include "utopia_petsc_Matrix.hpp"

#include "utopia_AppBase.hpp"

#ifdef USE_SIMD_ASSEMBLY
// #define USE_SIMD_LINEAR_ELASTICITY_FE
#endif

namespace utopia {

    template <class FunctionSpace>
    class LinearElasticityFE final : public Operator<typename FunctionSpace::Vector>, public Configurable {
    public:
        // using Mesh             = typename FunctionSpace::Mesh;
        using Elem = typename FunctionSpace::Elem;
        using Dev = typename FunctionSpace::Device;
        using Vector = typename FunctionSpace::Vector;
        using Matrix = typename FunctionSpace::Matrix;
        using Comm = typename FunctionSpace::Comm;

        static const int Dim = Elem::Dim;
        static const int NFunctions = Elem::NFunctions;

        using Point = typename FunctionSpace::Point;
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;
        using ElementMatrix = utopia::StaticMatrix<Scalar, NFunctions, NFunctions>;
        using ElementVector = utopia::StaticVector<Scalar, NFunctions>;

#ifdef USE_SIMD_LINEAR_ELASTICITY_FE
        using SIMDType = Vc::Vector<Scalar>;
        // using SIMDType = Scalar;
        using Quadrature = simd::Quadrature<SIMDType, Elem::Dim>;
        using GradValue = typename simd::FETraits<Elem, SIMDType>::GradValue;
#else
        using Quadrature = utopia::Quadrature<Elem, 2 * (Elem::Order - 1)>;
        using GradValue = typename Elem::GradValue;
#endif  // USE_SIMD_LINEAR_ELASTICITY_FE

        using Coefficient = utopia::Coefficient<FunctionSpace>;
        using LEKernel = utopia::LinearElasticityKernel<Scalar>;

        void read(Input &in) override {
            in.get("debug", debug_);
            in.get("mu", mu_);
            in.get("lambda", lambda_);
        }

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
            UTOPIA_TRACE_REGION_BEGIN("LinearElasticityFE::apply(...)");

            // const Comm &comm = space_.comm();

            if (y.empty()) {
                space_.create_vector(y);
            } else {
                y.set(0.0);
            }

            auto &&space_view = space_.view_device();

            x_coeff_->update(x);

            // LinearElasticity<FunctionSpace, Quadrature> elast(space_, quadrature, mu_, lambda_);
            Strain<FunctionSpace, Quadrature> strain(space_, quadrature_);

            {
                auto x_view = x_coeff_->view_device();
                auto y_view = space_.assembly_view_device(y);
                auto strain_view = strain.view_device();
                auto dx_view = differential_->view_device();

                Dev::parallel_for(
                    space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        Elem e;
                        ElementVector coeff, el_vec;
                        el_vec.set(0.0);

                        space_view.elem(i, e);
                        x_view.get(e, coeff);

                        auto &&strain = strain_view.make(e);
                        auto dx = dx_view.make(e);

                        const SizeType n_qp = strain.n_points();
                        const SizeType n_fun = strain.n_functions();

                        assert(n_qp > 0);

                        // taking advantage of symmetry
                        for (SizeType qp = 0; qp < n_qp; ++qp) {
                            const auto dx_qp = dx(qp);

                            for (SizeType j = 0; j < n_fun; ++j) {
                                auto &&strain_test = strain(j, qp);

                                const auto v_jj = simd::integrate(
                                    LEKernel::strain_apply(mu_, lambda_, strain_test, strain_test, dx_qp) * coeff(j));

                                assert(v_jj == v_jj);

                                el_vec(j) += v_jj;

                                for (SizeType l = j + 1; l < n_fun; ++l) {
                                    auto &&strain_trial = strain(l, qp);

                                    const auto v = simd::integrate(
                                        LEKernel::strain_apply(mu_, lambda_, strain_trial, strain_test, dx_qp));

                                    assert(v == v);

                                    el_vec(l) += v * coeff(j);
                                    el_vec(j) += v * coeff(l);
                                }
                            }
                        }

                        space_view.add_vector(e, el_vec, y_view);
                    });
            }

            space_.copy_at_constrained_dofs(x, y);
            assert(check_with_matrix(x, y));

            UTOPIA_TRACE_REGION_END("LinearElasticityFE::apply(...)");
            return true;
        }

        bool check_with_matrix(const Vector &x, const Vector &y) const {
            if (H_ptr_->empty()) {
                space_.create_matrix(*H_ptr_);
                hessian(x, *H_ptr_);
            }

            Vector mat_y = (*H_ptr_) * x;

            Scalar norm_diff = norm_infty(y - mat_y);

            // utopia::out() <<"norm_diff: " << norm_diff << std::endl;
            assert(norm_diff < 1e-10);
            return norm_diff < 1e-10;
        }

        bool hessian(const Vector &, Matrix &H) const {
            LinearElasticity<FunctionSpace, Quadrature> elast(space_, quadrature_, mu_, lambda_);
            MassMatrix<FunctionSpace, Quadrature> mass_matrix(space_, quadrature_);

            {
                auto &&space_view = space_.view_device();
                auto mat_view = space_.assembly_view_device(H);
                auto elast_view = elast.view_device();

                if (debug_) {
                    disp("elast");
                    elast_view.describe();
                }

                Dev::parallel_for(
                    space_.element_range(), UTOPIA_LAMBDA(const SizeType &i) {
                        Elem e;

                        // FIXME this is too big for GPU stack memory for hexas
                        ElementMatrix el_mat;
                        space_view.elem(i, e);

                        // Assemble local elast
                        el_mat.set(0.0);
                        elast_view.assemble(i, e, el_mat);
                        space_view.add_matrix(e, el_mat, mat_view);
                    });
            }

            space_.apply_constraints(H);
            return true;
        }

        LinearElasticityFE(FunctionSpace &space, const Scalar &mu = 1.0, const Scalar &lambda = 1.0)
            : space_(space),
              mu_(mu),
              lambda_(lambda),
              debug_(false),
              x_coeff_(nullptr),
              grad_(nullptr),
              differential_(nullptr) {
            H_ptr_ = utopia::make_unique<Matrix>();
            init();
        }

    private:
        FunctionSpace space_;
        Scalar mu_, lambda_;
        bool debug_;

        Quadrature quadrature_;
        std::unique_ptr<Coefficient> x_coeff_;
        std::unique_ptr<PhysicalGradient<FunctionSpace, Quadrature>> grad_;
        std::unique_ptr<Differential<FunctionSpace, Quadrature>> differential_;

        std::unique_ptr<Matrix> H_ptr_;

        void init() {
            x_coeff_ = utopia::make_unique<Coefficient>(space_);

#ifdef USE_SIMD_LINEAR_ELASTICITY_FE
            utopia::out() << "Using SIMD based LinearElasticityFE\n";
            simd::QuadratureDB<Elem, SIMDType>::get(2 * (Elem::Order - 1), quadrature_);
            utopia::out() << "simd lanes: " << SIMDType::Size << ", qp: " << quadrature_.n_points() << '\n';
#endif  // USE_SIMD_LINEAR_ELASTICITY_FE

            grad_ = utopia::make_unique<PhysicalGradient<FunctionSpace, Quadrature>>(space_, quadrature_);
            differential_ = utopia::make_unique<Differential<FunctionSpace, Quadrature>>(space_, quadrature_);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_LINEAR_ELASTICITY_FE_HPP

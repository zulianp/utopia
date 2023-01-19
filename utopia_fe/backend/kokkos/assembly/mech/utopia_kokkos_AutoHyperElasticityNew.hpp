#ifndef UTOPIA_KOKKOS_AUTOHYPERELASTICITY_NEW_HPP
#define UTOPIA_KOKKOS_AUTOHYPERELASTICITY_NEW_HPP

#include "utopia_fe_base.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Field.hpp"
#include "utopia_kokkos_Gradient.hpp"
#include "utopia_kokkos_Material.hpp"
#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE_, class Material>
        class AutoHyperElasticityNew : public utopia::Material<FunctionSpace, FE_> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using Super = utopia::Material<FunctionSpace, FE_>;
            using Gradient = utopia::kokkos::Gradient<FE>;
            using DynRankView = typename Gradient::DynRankView;

            static constexpr int Dim = Material::Dim;
            using Params = typename Material::Params;

            AutoHyperElasticityNew(Params op = Params()) : Super(), material_(std::move(op)) {}

            int order() const override { return order_; }

            inline int n_vars() const override { return Dim; }

            inline bool has_hessian() const override { return true; }
            inline bool has_gradient() const override { return true; }
            inline bool has_value() const override { return true; }
            inline bool is_operator() const override { return false; }
            inline bool is_linear() const override { return false; }

            inline std::string name() const override {
                return std::string("AutoHyperElasticityNew<") + Material::class_name() + ">";
            }

            std::shared_ptr<Gradient> compute_deformation_gradient() const {
                auto displacement = this->assembler()->current_solution();

                assert(displacement);
                assert(displacement->is_coefficient());

                if (!displacement->is_coefficient()) {
                    Utopia::Abort("compute_deformation_gradient: displacement must me in coefficient form!");
                }

                auto deformation_gradient = std::make_shared<Gradient>(this->assembler()->fe_ptr());

                deformation_gradient->init(*displacement);
                deformation_gradient->add_identity();
                assert(deformation_gradient->check_dets_are_positive());
                return deformation_gradient;
            }

            bool apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) override { return false; }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Value
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            template <class Measure>
            class ValueKernel {
            public:
                UTOPIA_INLINE_FUNCTION ValueKernel(const Measure &measure,
                                                   const Material &material,
                                                   const DynRankView &deformation_gradient,
                                                   int n_qp)
                    : measure(measure), material(material), deformation_gradient(deformation_gradient), n_qp(n_qp) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell) const {
                    StaticMatrix<Scalar, Dim, Dim> F;

                    Scalar ret = 0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        const Scalar dx = measure(cell, qp);

                        for (int d1 = 0; d1 < Dim; ++d1) {
                            for (int d2 = 0; d2 < Dim; ++d2) {
                                F(d1, d2) = deformation_gradient(cell, qp, d1 * Dim + d2);
                            }
                        }

                        Scalar val = 0;
                        material.value(&F.raw_type()[0], dx, val);
                        ret += val;
                    }

                    return ret;
                }

                Measure measure;
                Material material;
                DynRankView deformation_gradient;
                int n_qp;
            };

            template <class Measure>
            inline ValueKernel<Measure> value_kernel(const Measure &measure,
                                                     const Material &material,
                                                     const DynRankView &deformation_gradient,
                                                     int n_qp) {
                return ValueKernel<Measure>(measure, material, deformation_gradient, n_qp);
            }

            bool value_assemble(AssemblyMode mode) override {
                UTOPIA_TRACE_REGION_BEGIN(name() + "::value_assemble");

                auto &&assembler = this->assembler();
                assert(assembler);
                auto &&fe = assembler->fe();

                auto deformation_gradient = this->compute_deformation_gradient();

                auto k = value_kernel(fe.measure(), material_, deformation_gradient->data(), fe.n_quad_points());

                assembler->assemble_scalar_e(name() + "::value", mode, k);

                UTOPIA_TRACE_REGION_END(name() + "::value_assemble");
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            // Hessian
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            template <class Grad, class Measure>
            class HessianKernel : public TestTrialOp {
            public:
                static const int NComponentsTest = Dim;
                static const int NComponentsTrial = Dim;

                UTOPIA_INLINE_FUNCTION HessianKernel(const Grad &grad,
                                                     const Measure &measure,
                                                     const Material &material,
                                                     const DynRankView &deformation_gradient,
                                                     int n_qp)
                    : grad(grad),
                      measure(measure),
                      material(material),
                      deformation_gradient(deformation_gradient),
                      n_qp(n_qp) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                       const int i,
                                                       const int j,
                                                       StaticMatrix<Scalar, Dim, Dim> &block) const {
                    StaticVector<Scalar, Dim> grad_test, grad_trial;
                    StaticMatrix<Scalar, Dim, Dim> F;

                    block.set(0.0);

                    for (int qp = 0; qp < n_qp; ++qp) {
                        for (int d1 = 0; d1 < Dim; ++d1) {
                            for (int d2 = 0; d2 < Dim; ++d2) {
                                F(d1, d2) = deformation_gradient(cell, qp, d1 * Dim + d2);
                            }
                        }

                        Scalar dx = measure(cell, qp);

                        for (int d = 0; d < Dim; ++d) {
                            grad_test[d] = grad(cell, i, qp, d);
                        }

                        for (int d = 0; d < Dim; ++d) {
                            grad_trial[d] = grad(cell, j, qp, d);
                        }

                        material.hessian(&F.raw_type()[0], &grad_test[0], &grad_trial[0], dx, &block.raw_type()[0]);
                    }
                }

                Grad grad;
                Measure measure;
                Material material;
                DynRankView deformation_gradient;
                int n_qp;
            };

            template <class Grad, class Measure>
            inline HessianKernel<Grad, Measure> hessian_kernel(const Grad &grad,
                                                               const Measure &measure,
                                                               const Material &material,
                                                               const DynRankView &deformation_gradient,
                                                               int n_qp) {
                return HessianKernel<Grad, Measure>(grad, measure, material, deformation_gradient, n_qp);
            }

            bool hessian_assemble(AssemblyMode mode) override {
                UTOPIA_TRACE_REGION_BEGIN(name() + "::hessian_assemble");

                auto &&assembler = this->assembler();
                assert(assembler);
                auto &&fe = assembler->fe();

                auto deformation_gradient = this->compute_deformation_gradient();

                auto k = hessian_kernel(
                    fe.grad(), fe.measure(), material_, deformation_gradient->data(), fe.n_quad_points());

                assembler->assemble_matrix_eij_block(name() + "::hessian", mode, k);

                // std::cout << "------------------------------------------------\n";
                // std::cout << "------------------------------------------------\n";
                // const SizeType e0 = assembler->matrix_accumulator()->data().extent(0);
                // const SizeType e1 = assembler->matrix_accumulator()->data().extent(1);
                // const SizeType e2 = assembler->matrix_accumulator()->data().extent(2);
                // const SizeType e3 = assembler->matrix_accumulator()->data().extent(2);

                // std::cout << "H: " << e0 << " " << e1 << " " << e2 << " " << e3 << "\n";

                // assembler->matrix_accumulator()->describe(std::cout);

                // std::cout << "------------------------------------------------\n";

                UTOPIA_TRACE_REGION_END(name() + "::hessian_assemble");
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            /// Gradient
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            template <class Grad, class Measure>
            class GradientKernel : public TestOp {
            public:
                static const int NComponentsTest = Dim;

                UTOPIA_INLINE_FUNCTION GradientKernel(const Grad &grad,
                                                      const Measure &measure,
                                                      const Material &material,
                                                      const DynRankView &deformation_gradient,
                                                      int n_qp)
                    : grad(grad),
                      measure(measure),
                      material(material),
                      deformation_gradient(deformation_gradient),
                      n_qp(n_qp) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell,
                                                       const int i,
                                                       StaticVector<Scalar, Dim> &block) const {
                    StaticVector<Scalar, Dim> grad_test;
                    StaticMatrix<Scalar, Dim, Dim> F;

                    block.set(0.0);

                    for (int qp = 0; qp < n_qp; ++qp) {
                        for (int d1 = 0; d1 < Dim; ++d1) {
                            for (int d2 = 0; d2 < Dim; ++d2) {
                                F(d1, d2) = deformation_gradient(cell, qp, d1 * Dim + d2);
                            }
                        }

                        Scalar dx = measure(cell, qp);

                        for (int d = 0; d < Dim; ++d) {
                            grad_test[d] = grad(cell, i, qp, d);
                        }

                        material.gradient(&F.raw_type()[0], &grad_test[0], dx, &block.raw_type()[0]);
                    }
                }

                Grad grad;
                Measure measure;
                Material material;
                DynRankView deformation_gradient;
                int n_qp;
            };

            template <class Grad, class Measure>
            inline GradientKernel<Grad, Measure> gradient_kernel(const Grad &grad,
                                                                 const Measure &measure,
                                                                 const Material &material,
                                                                 const DynRankView &deformation_gradient,
                                                                 int n_qp) {
                return GradientKernel<Grad, Measure>(grad, measure, material, deformation_gradient, n_qp);
            }

            bool gradient_assemble(AssemblyMode mode) override {
                UTOPIA_TRACE_REGION_BEGIN(name() + "::gradient_assemble");

                auto &&assembler = this->assembler();
                assert(assembler);
                auto &&fe = assembler->fe();

                auto deformation_gradient = this->compute_deformation_gradient();

                // std::cout << "------------------------------------------------\n";
                // std::cout << "------------------------------------------------\n";
                // const SizeType e0 = deformation_gradient->data().extent(0);
                // const SizeType e1 = deformation_gradient->data().extent(1);
                // const SizeType e2 = deformation_gradient->data().extent(2);

                // std::cout << "F: " << e0 << " " << e1 << " " << e2 << "\n";

                // deformation_gradient->describe(std::cout);

                // std::cout << "------------------------------------------------\n";

                auto k = gradient_kernel(
                    fe.grad(), fe.measure(), material_, deformation_gradient->data(), fe.n_quad_points());

                assembler->assemble_vector_ei_block(name() + "::gradient", mode, k);

                // std::cout << "------------------------------------------------\n";
                // std::cout << "------------------------------------------------\n";
                // const SizeType e0 = assembler->vector_accumulator()->data().extent(0);
                // const SizeType e1 = assembler->vector_accumulator()->data().extent(1);
                // const SizeType e2 = assembler->vector_accumulator()->data().extent(2);

                // std::cout << "G: " << e0 << " " << e1 << " " << e2 << "\n";

                // assembler->vector_accumulator()->describe(std::cout);

                // std::cout << "------------------------------------------------\n";

                UTOPIA_TRACE_REGION_END(name() + "::gradient_assemble");
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            // NVCC_PRIVATE :
            Material material_;
            int order_{6};
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_AUTOHYPERELASTICITY_NEW_HPP

#ifndef UTOPIA_KOKKOS_MIXED_HYPER_ELASTIC_MATERIAL_HPP
#define UTOPIA_KOKKOS_MIXED_HYPER_ELASTIC_MATERIAL_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Gradient.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    namespace kokkos {

        template <class FE_, class Material>
        class MixedHyperElasticity : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FE_, DefaultView<typename FE_::Scalar>>;

            static constexpr int Dim = Material::Dim;
            static constexpr int block_size = Dim + 1;

            using Params = typename Material::Params;

            MixedHyperElasticity(const std::shared_ptr<FE> &fe, Params op = Params())
                : Super(fe), material_(std::move(op)) {
                assert(!fe || fe->spatial_dimension() == Dim);
            }

            inline int n_vars() const override { return block_size; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            inline bool is_linear() const override { return false; }
            inline bool is_operator() const override { return false; }

            inline std::string name() const override { return "MixedHyperElasticity"; }

            virtual bool update(const std::shared_ptr<Field<FE>> &field) override {
                if (!Super::update(field)) {
                    return false;
                }

                assert(field);
                assert(field->is_coefficient());

                if (!field->is_coefficient()) {
                    Utopia::Abort("Assemble<MixedHyperElasticity>::update, field must me in coefficient form!");
                }

                if (!deformation_gradient_) {
                    // Initialize gradient
                    deformation_gradient_ = std::make_shared<Gradient<FE>>(this->fe_ptr());
                }

                field_ = field;

                deformation_gradient_->init(*field, 0, Dim);
                deformation_gradient_->add_identity();
                assert(deformation_gradient_->check_dets_are_positive());
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<MixedHyperElasticity>::assemble_matrix");

                this->ensure_matrix_accumulator();

                auto &fe = this->fe();
                auto data = this->matrix_data();

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();

                    auto displacement_and_pressure = field_->interpolate();

                    auto material = material_;

                    int n_quad_points = this->fe().n_quad_points();
                    int n_shape_functions = this->fe().n_shape_functions();

                    auto grad = this->fe().grad();
                    auto fun = this->fe().fun();

                    auto measure = this->fe().measure();

                    this->loop_cell(
                        "Hessian", UTOPIA_LAMBDA(int cell) {
                            StaticVector<Scalar, Dim> grad_test, grad_trial;
                            StaticMatrix<Scalar, Dim, Dim> F_qp;
                            StaticMatrix<Scalar, Dim+1, Dim+1> stress;

                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                // Copy deformation gradient at qp
                                for (int d1 = 0; d1 < Dim; ++d1) {
                                    for (int d2 = 0; d2 < Dim; ++d2) {
                                        F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
                                    }
                                }

                                // Interpolate pressure
                                Scalar p = displacement_and_pressure(cell, qp, Dim);

                                Scalar dx = measure(cell, qp);

                                for (int i = 0; i < n_shape_functions; ++i) {
                                    // Copy deformation gradient at qp
                                    for (int d1 = 0; d1 < Dim; ++d1) {
                                        grad_test[d1] = grad(cell, i, qp, d1);
                                    }

                                    for (int j = 0; j < n_shape_functions; ++j) {
                                        for (int d1 = 0; d1 < Dim; ++d1) {
                                            grad_trial[d1] = grad(cell, j, qp, d1);
                                        }
                                        // Assemble
                                        stress.set(0.);

                                        material.hessian(&F_qp.raw_type()[0],
                                                         p,
                                                         &grad_test[0],
                                                         &grad_trial[0],
                                                         fun(i, qp),
                                                         fun(j, qp),
                                                         dx,
                                                         &stress.raw_type()[0]);

                                        for (int d1 = 0; d1 < block_size; ++d1) {
                                            auto dof_i = i * block_size + d1;
                                            for (int d2 = 0; d2 < block_size; ++d2) {
                                                auto dof_j = j * block_size + d2;

                                                assert(stress(d1, d2) == stress(d1, d2));
                                                data(cell, dof_i, dof_j) += stress(d1, d2);
                                            }
                                        }
                                    }
                                }
                            }
                        });
                }

                UTOPIA_TRACE_REGION_END("Assemble<MixedHyperElasticity>::assemble_matrix");
                return true;
            }

            bool assemble_vector() override {
                UTOPIA_TRACE_REGION_BEGIN("Assemble<MixedHyperElasticity>::assemble_vector");

                this->ensure_vector_accumulator();

                auto &fe = this->fe();
                auto data = this->vector_data();

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();

                    auto displacement_and_pressure = field_->interpolate();

                    auto material = material_;

                    int n_quad_points = this->fe().n_quad_points();
                    int n_shape_functions = this->fe().n_shape_functions();

                    auto grad = this->fe().grad();
                    auto fun = this->fe().fun();
                    auto measure = this->fe().measure();

                    this->loop_cell(
                        "Gradient", UTOPIA_LAMBDA(int cell) {
                            StaticVector<Scalar, Dim> grad_test;
                            StaticMatrix<Scalar, Dim, Dim> F_qp;
                            StaticVector<Scalar, Dim+1> stress;

                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                // Copy deformation gradient at qp
                                for (int d1 = 0; d1 < Dim; ++d1) {
                                    for (int d2 = 0; d2 < Dim; ++d2) {
                                        F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
                                    }
                                }

                                // Interpolate pressure
                                Scalar p = displacement_and_pressure(cell, qp, Dim);

                                Scalar dx = measure(cell, qp);

                                for (int i = 0; i < n_shape_functions; ++i) {
                                    // Copy deformation gradient at qp
                                    for (int d1 = 0; d1 < Dim; ++d1) {
                                        grad_test[d1] = grad(cell, i, qp, d1);
                                    }

                                    // Assemble
                                    stress.set(0.);
                                    material.gradient(
                                        &F_qp.raw_type()[0], p, &grad_test[0], fun(i, qp), dx, &stress.raw_type()[0]);

                                    for (int d1 = 0; d1 < block_size; ++d1) {
                                        auto dof_i = i * block_size + d1;

                                        assert(stress[d1] == stress[d1]);
                                        data(cell, dof_i) += stress[d1];
                                    }
                                }
                            }
                        });
                }

                UTOPIA_TRACE_REGION_END("Assemble<MixedHyperElasticity>::assemble_vector");
                return true;
            }

            // NVCC_PRIVATE :
            Material material_;
            std::shared_ptr<Gradient<FE>> deformation_gradient_;
            std::shared_ptr<Field<FE>> field_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_MIXED_HYPER_ELASTIC_MATERIAL_HPP

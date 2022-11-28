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

            static constexpr int Dim = Material::Dim;
            using Params = typename Material::Params;

            AutoHyperElasticityNew(Params op = Params()) : Super(), material_(std::move(op)) {}

            inline int n_vars() const override { return Dim; }

            bool has_hessian() const override { return true; }
            bool has_gradient() const override { return true; }
            bool has_value() const override { return true; }
            bool is_operator() const override { return true; }

            inline std::string name() const override {
                return std::string("AutoHyperElasticityNew<") + Material::class_name() + ">";
            }

            bool hessian_assemble(AssemblyMode mode) override { return false; }
            bool gradient_assemble(AssemblyMode) override { return false; }
            bool value_assemble(AssemblyMode mode) override { return false; }
            bool apply_assemble(utopia::kokkos::Field<FE> &field, AssemblyMode mode) override { return false; }

            // virtual bool update(const std::shared_ptr<Field<FE>> &displacement) override {
            //     if (!Super::update(displacement)) {
            //         return false;
            //     }

            //     assert(displacement);
            //     assert(displacement->is_coefficient());

            //     if (!displacement->is_coefficient()) {
            //         Utopia::Abort(name() + ":update, displacement must me in coefficient form!");
            //     }

            //     if (!deformation_gradient_) {
            //         // Initialize gradient
            //         deformation_gradient_ = std::make_shared<Gradient<FE>>(this->fe_ptr());
            //     }

            //     deformation_gradient_->init(*displacement);
            //     deformation_gradient_->add_identity();
            //     assert(deformation_gradient_->check_dets_are_positive());
            //     return true;
            // }

            // bool assemble_matrix() override {
            //     UTOPIA_TRACE_REGION_BEGIN(name() + "::assemble_matrix");

            //     this->ensure_matrix_accumulator();

            //     auto &fe = this->fe();
            //     auto data = this->matrix_data();

            //     {
            //         assert(deformation_gradient_);
            //         auto F = deformation_gradient_->data();
            //         auto material = material_;

            //         int n_quad_points = this->fe().n_quad_points();
            //         int n_shape_functions = this->fe().n_shape_functions();

            //         auto grad = this->fe().grad();
            //         auto measure = this->fe().measure();

            //         this->loop_cell(
            //             "Hessian", UTOPIA_LAMBDA(int cell) {
            //                 StaticVector<Scalar, Dim> grad_test, grad_trial;
            //                 StaticMatrix<Scalar, Dim, Dim> F_qp, stress;

            //                 for (int qp = 0; qp < n_quad_points; ++qp) {
            //                     // Copy deformation gradient at qp
            //                     for (int d1 = 0; d1 < Dim; ++d1) {
            //                         for (int d2 = 0; d2 < Dim; ++d2) {
            //                             F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
            //                         }
            //                     }

            //                     Scalar dx = measure(cell, qp);

            //                     for (int i = 0; i < n_shape_functions; ++i) {
            //                         // Copy deformation gradient at qp
            //                         for (int d1 = 0; d1 < Dim; ++d1) {
            //                             grad_test[d1] = grad(cell, i, qp, d1);
            //                         }

            //                         for (int j = 0; j < n_shape_functions; ++j) {
            //                             for (int d1 = 0; d1 < Dim; ++d1) {
            //                                 grad_trial[d1] = grad(cell, j, qp, d1);
            //                             }
            //                             // Assemble
            //                             stress.set(0.);

            //                             material.hessian(&F_qp.raw_type()[0],
            //                                              &grad_test[0],
            //                                              &grad_trial[0],
            //                                              dx,
            //                                              &stress.raw_type()[0]);

            //                             for (int d1 = 0; d1 < Dim; ++d1) {
            //                                 auto dof_i = i * Dim + d1;
            //                                 for (int d2 = 0; d2 < Dim; ++d2) {
            //                                     auto dof_j = j * Dim + d2;

            //                                     assert(stress(d1, d2) == stress(d1, d2));
            //                                     data(cell, dof_i, dof_j) += stress(d1, d2);
            //                                 }
            //                             }
            //                         }
            //                     }
            //                 }
            //             });
            //     }

            //     UTOPIA_TRACE_REGION_END(name() + "::assemble_matrix");
            //     return true;
            // }

            // bool assemble_vector() override {
            //     UTOPIA_TRACE_REGION_BEGIN(name() + "::assemble_vector");

            //     this->ensure_vector_accumulator();

            //     auto &fe = this->fe();
            //     auto data = this->vector_data();

            //     {
            //         assert(deformation_gradient_);
            //         auto F = deformation_gradient_->data();

            //         auto material = material_;

            //         int n_quad_points = this->fe().n_quad_points();
            //         int n_shape_functions = this->fe().n_shape_functions();

            //         auto grad = this->fe().grad();
            //         auto measure = this->fe().measure();

            //         this->loop_cell(
            //             "Gradient", UTOPIA_LAMBDA(int cell) {
            //                 StaticVector<Scalar, Dim> stress, grad_test;
            //                 StaticMatrix<Scalar, Dim, Dim> F_qp;

            //                 for (int qp = 0; qp < n_quad_points; ++qp) {
            //                     // Copy deformation gradient at qp
            //                     for (int d1 = 0; d1 < Dim; ++d1) {
            //                         for (int d2 = 0; d2 < Dim; ++d2) {
            //                             F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
            //                         }
            //                     }

            //                     Scalar dx = measure(cell, qp);

            //                     for (int i = 0; i < n_shape_functions; ++i) {
            //                         // Copy deformation gradient at qp
            //                         for (int d1 = 0; d1 < Dim; ++d1) {
            //                             grad_test[d1] = grad(cell, i, qp, d1);
            //                         }

            //                         // Assemble
            //                         stress.set(0.);
            //                         material.gradient(&F_qp.raw_type()[0], &grad_test[0], dx, &stress.raw_type()[0]);

            //                         for (int d1 = 0; d1 < Dim; ++d1) {
            //                             auto dof_i = i * Dim + d1;

            //                             assert(stress[d1] == stress[d1]);
            //                             data(cell, dof_i) += stress[d1];
            //                         }
            //                     }
            //                 }
            //             });
            //     }

            //     UTOPIA_TRACE_REGION_END(name() + "::assemble_vector");
            //     return true;
            // }

            // NVCC_PRIVATE :
            Material material_;
            std::shared_ptr<Gradient<FE>> deformation_gradient_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_AUTOHYPERELASTICITY_NEW_HPP

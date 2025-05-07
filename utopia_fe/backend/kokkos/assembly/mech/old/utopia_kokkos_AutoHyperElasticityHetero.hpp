#ifndef UTOPIA_KOKKOS_AUTOHYPERELASTICITY_HETERO_HPP
#define UTOPIA_KOKKOS_AUTOHYPERELASTICITY_HETERO_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Gradient.hpp"

#include "utopia_StressStrainParameters.hpp"
#include "utopia_Views.hpp"

namespace utopia {
    namespace kokkos {

        template <class FunctionSpace, class FE_, class Material>
        class AutoHyperElasticityHetero : public FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;

            using ExecutionSpace = typename FE::ExecutionSpace;
            using Super = utopia::kokkos::FEAssembler<FunctionSpace, FE_, DefaultView<typename FE_::Scalar>>;

            static constexpr int Dim = Material::Dim;

            class Params : public Configurable {
            public:
                using MaterialParams = typename Material::Params;

                static const int MAX_SUBDOMAINS = 10;
                MaterialParams default_params;

                int n_params{0};
                int blocks[MAX_SUBDOMAINS];
                MaterialParams params[MAX_SUBDOMAINS];

                void read(Input &in) override {
                    bool verbose = false;
                    in.get("verbose", verbose);
                    default_params.read(in);

                    in.get("blocks", [this](Input &array_node) {
                        array_node.get_all([this](Input &node) {
                            int block;
                            node.require("block", block);

                            bool verbose = false;
                            node.get("verbose", verbose);
                            if (verbose && !mpi_world_rank()) {
                                utopia::out() << "block: " << block << "\n";
                            }

                            MaterialParams p;
                            p.read(node);

                            assert(n_params < MAX_SUBDOMAINS);
                            blocks[n_params] = block;
                            params[n_params] = p;
                            n_params++;
                        });
                    });
                }
            };

            class BlockParams {
            public:
                using MaterialParams = typename Material::Params;

                BlockParams(const typename FE::IntView &tags, const Params &params) : tags(tags), params(params) {}
                UTOPIA_FUNCTION MaterialParams operator()(const int cell) const {
                    auto t = tags(cell);

                    MaterialParams ret = params.default_params;

                    for (int i = 0; i < params.n_params; i++) {
                        if (t == params.blocks[i]) {
                            ret = params.params[i];
                            break;
                        }
                    }

                    return ret;
                }

                typename FE::IntView tags;
                Params params;
            };

            AutoHyperElasticityHetero(const std::shared_ptr<FE> &fe, Params params = Params())
                : Super(fe), mapped_params_(params) {
                assert(!fe || fe->spatial_dimension() == Dim);
            }

            inline int n_vars() const override { return this->fe().spatial_dimension(); }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            inline bool is_linear() const override { return false; }
            inline bool is_operator() const override { return false; }

            inline std::string name() const override {
                return std::string("AutoHyperElasticityHetero<") + Material::class_name() + ">";
            }

            virtual bool update(const std::shared_ptr<Field<FE>> &displacement) override {
                if (!Super::update(displacement)) {
                    return false;
                }

                assert(displacement);
                assert(displacement->is_coefficient());

                if (!displacement->is_coefficient()) {
                    Utopia::Abort(name() + ":update, displacement must me in coefficient form!");
                }

                if (!deformation_gradient_) {
                    // Initialize gradient
                    deformation_gradient_ = std::make_shared<Gradient<FE>>(this->fe_ptr());
                }

                deformation_gradient_->init(*displacement);
                deformation_gradient_->add_identity();
                assert(deformation_gradient_->check_dets_are_positive());

                // std::cout << "------------------------------------------------\n";
                // std::cout << "------------------------------------------------\n";
                // const SizeType e0 = deformation_gradient_->data().extent(0);
                // const SizeType e1 = deformation_gradient_->data().extent(1);
                // const SizeType e2 = deformation_gradient_->data().extent(2);

                // std::cout << "F: " << e0 << " " << e1 << " " << e2 << "\n";

                // deformation_gradient_->describe(std::cout);

                // std::cout << "------------------------------------------------\n";
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN(name() + "::assemble_matrix");

                this->ensure_matrix_accumulator();

                auto &fe = this->fe();
                auto data = this->matrix_data();

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();
                    auto material = material_;

                    int n_quad_points = fe.n_quad_points();
                    int n_shape_functions = fe.n_shape_functions();

                    auto grad = fe.grad();
                    auto measure = fe.measure();

                    assert(fe.has_element_tags());
                    if (!fe.has_element_tags()) {
                        Utopia::Abort("Material requires element tags!");
                    }

                    BlockParams block_params(fe.element_tags(), mapped_params_);

                    this->loop_cell(
                        "Hessian", UTOPIA_LAMBDA(int cell) {
                            StaticVector<Scalar, Dim> grad_test, grad_trial;
                            StaticMatrix<Scalar, Dim, Dim> F_qp, stress;

                            auto params = block_params(cell);

                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                // Copy deformation gradient at qp
                                for (int d1 = 0; d1 < Dim; ++d1) {
                                    for (int d2 = 0; d2 < Dim; ++d2) {
                                        F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
                                    }
                                }

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

                                        material.hessian(params,
                                                         &F_qp.raw_type()[0],
                                                         &grad_test[0],
                                                         &grad_trial[0],
                                                         dx,
                                                         &stress.raw_type()[0]);

                                        for (int d1 = 0; d1 < Dim; ++d1) {
                                            auto dof_i = i * Dim + d1;
                                            for (int d2 = 0; d2 < Dim; ++d2) {
                                                auto dof_j = j * Dim + d2;

                                                assert(stress(d1, d2) == stress(d1, d2));
                                                data(cell, dof_i, dof_j) += stress(d1, d2);
                                            }
                                        }
                                    }
                                }
                            }
                        });

                    // std::cout << "------------------------------------------------\n";
                    // std::cout << "------------------------------------------------\n";
                    // const SizeType e0 = this->matrix_accumulator()->data().extent(0);
                    // const SizeType e1 = this->matrix_accumulator()->data().extent(1);
                    // const SizeType e2 = this->matrix_accumulator()->data().extent(2);
                    // const SizeType e3 = this->matrix_accumulator()->data().extent(2);

                    // std::cout << "H: " << e0 << " " << e1 << " " << e2 << " " << e3 << "\n";

                    // this->matrix_accumulator()->describe(std::cout);

                    // std::cout << "------------------------------------------------\n";
                }

                UTOPIA_TRACE_REGION_END(name() + "::assemble_matrix");
                return true;
            }

            bool assemble_vector() override {
                UTOPIA_TRACE_REGION_BEGIN(name() + "::assemble_vector");

                this->ensure_vector_accumulator();

                auto &fe = this->fe();
                auto data = this->vector_data();

                {
                    assert(deformation_gradient_);
                    auto F = deformation_gradient_->data();

                    auto material = material_;

                    int n_quad_points = this->fe().n_quad_points();
                    int n_shape_functions = this->fe().n_shape_functions();

                    auto grad = this->fe().grad();
                    auto measure = this->fe().measure();

                    assert(fe.has_element_tags());
                    if (!fe.has_element_tags()) {
                        Utopia::Abort("Material requires element tags!");
                    }

                    BlockParams block_params(fe.element_tags(), mapped_params_);

                    this->loop_cell(
                        "Gradient", UTOPIA_LAMBDA(int cell) {
                            StaticVector<Scalar, Dim> stress, grad_test;
                            StaticMatrix<Scalar, Dim, Dim> F_qp;

                            auto params = block_params(cell);

                            for (int qp = 0; qp < n_quad_points; ++qp) {
                                // Copy deformation gradient at qp
                                for (int d1 = 0; d1 < Dim; ++d1) {
                                    for (int d2 = 0; d2 < Dim; ++d2) {
                                        F_qp(d1, d2) = F(cell, qp, d1 * Dim + d2);
                                    }
                                }

                                Scalar dx = measure(cell, qp);

                                for (int i = 0; i < n_shape_functions; ++i) {
                                    // Copy deformation gradient at qp
                                    for (int d1 = 0; d1 < Dim; ++d1) {
                                        grad_test[d1] = grad(cell, i, qp, d1);
                                    }

                                    // Assemble
                                    stress.set(0.);
                                    material.gradient(
                                        params, &F_qp.raw_type()[0], &grad_test[0], dx, &stress.raw_type()[0]);

                                    for (int d1 = 0; d1 < Dim; ++d1) {
                                        auto dof_i = i * Dim + d1;

                                        assert(stress[d1] == stress[d1]);
                                        data(cell, dof_i) += stress[d1];
                                    }
                                }
                            }
                        });
                }

                // std::cout << "------------------------------------------------\n";
                // std::cout << "------------------------------------------------\n";
                // const SizeType e0 = this->vector_accumulator()->data().extent(0);
                // const SizeType e1 = this->vector_accumulator()->data().extent(1);
                // const SizeType e2 = this->vector_accumulator()->data().extent(2);

                // std::cout << "G: " << e0 << " " << e1 << " " << e2 << "\n";

                // this->vector_accumulator()->describe(std::cout);

                // std::cout << "------------------------------------------------\n";

                UTOPIA_TRACE_REGION_END(name() + "::assemble_vector");
                return true;
            }

            // NVCC_PRIVATE :
            Material material_;
            std::shared_ptr<Gradient<FE>> deformation_gradient_;
            Params mapped_params_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_AUTOHYPERELASTICITY_HETERO_HPP

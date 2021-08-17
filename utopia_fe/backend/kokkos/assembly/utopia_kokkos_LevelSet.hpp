#ifndef UTOPIA_KOKKOS_LEVELSET_HPP
#define UTOPIA_KOKKOS_LEVELSET_HPP

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Gradient.hpp"

namespace utopia {
    namespace kokkos {

        template <class FE_, class FirstLameParameter = typename FE_::Scalar, class ShearModulus = FirstLameParameter>
        class LevelSet : public FEAssembler<FE_, DefaultView<typename FE_::Scalar>> {
        public:
            using FE = FE_;
            using SizeType = typename FE::SizeType;
            using Scalar = typename FE::Scalar;
            using DynRankView = typename FE::DynRankView;
            using ExecutionSpace = typename FE::ExecutionSpace;

            using Super = utopia::kokkos::FEAssembler<FE_, DefaultView<typename FE_::Scalar>>;

            class Op {
            public:
                UTOPIA_INLINE_FUNCTION Op(const Field &gradient_field,
                                          const Grad &grad,
                                          const Fun &fun,
                                          const Measure &measure)
                    : gradient_field(gradient_field),
                      grad(grad),
                      fun(fun),
                      measure(measure),
                      n_qp(measure.extent(1)),
                      dim(grad.extent(3)) {
                    assert(gradient_field.extent(0) == measure.extent(0));
                    assert(gradient_field.extent(1) == measure.extent(1));
                    assert(gradient_field.extent(2) == grad.extent(3));
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i, const int &j) const {
                    Scalar integral = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar val = 0.0;
                        for (int dj = 0; dj < dim; ++dj) {
                            assert(gradient_field(cell, qp, dj) == gradient_field(cell, qp, dj));
                            assert(grad(cell, j, qp, dj) == grad(cell, j, qp, dj));

                            val += grad(cell, j, qp, dj) * gradient_field(cell, qp, dj);
                        }

                        assert(val == val);

                        auto dX = measure(cell, qp);
                        assert(dX == dX);

                        val *= fun(i, qp) * dX;
                        integral += val;
                    }

                    assert(integral == integral);
                    return 2 * integral;
                }

                Field gradient_field;
                Grad grad;
                Fun fun;
                Measure measure;
                const int n_qp;
                const int dim;
            };

            class OpStab {
            public:
                UTOPIA_INLINE_FUNCTION OpStab(const Field &gradient_field,
                                              const Grad &grad,
                                              const Fun &fun,
                                              const Measure &measure)
                    : gradient_field(gradient_field),
                      grad(grad),
                      fun(fun),
                      measure(measure),
                      n_qp(measure.extent(1)),
                      dim(grad.extent(3)) {
                    assert(gradient_field.extent(0) == measure.extent(0));
                    assert(gradient_field.extent(1) == measure.extent(1));
                    assert(gradient_field.extent(2) == grad.extent(3));
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i, const int &j) const {
                    Scalar integral = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar norm_grad = 0.0;

                        for (int dj = 0; dj < dim; ++dj) {
                            assert(gradient_field(cell, qp, dj) == gradient_field(cell, qp, dj));
                            norm_grad += gradient_field(cell, j, qp, dj) * gradient_field(cell, qp, dj);
                        }

                        norm_grad = device::sqrt(norm_grad);

                        assert(val == val);

                        auto dX = measure(cell, qp);
                        assert(dX == dX);

                        val *= fun(i, qp) * dX;
                        integral += val;
                    }

                    assert(integral == integral);
                    return 2 * integral;
                }

                Field gradient_field;
                Grad grad;
                Fun fun;
                Measure measure;
                const int n_qp;
                const int dim;
            };

            class GradientOp {
            public:
                UTOPIA_INLINE_FUNCTION Op(const Field &gradient_field, const Fun &fun, const Measure &measure)
                    : gradient_field(gradient_field),
                      fun(fun),
                      measure(measure),
                      n_qp(measure.extent(1)),
                      dim(grad.extent(3)) {
                    assert(gradient_field.extent(0) == measure.extent(0));
                    assert(gradient_field.extent(1) == measure.extent(1));
                    assert(gradient_field.extent(2) == grad.extent(3));
                }

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int &cell, const int &i) const {
                    Scalar integral = 0.0;
                    for (int qp = 0; qp < n_qp; ++qp) {
                        Scalar val = 0.0;
                        for (int dj = 0; dj < dim; ++dj) {
                            assert(gradient_field(cell, qp, dj) == gradient_field(cell, qp, dj));
                            val += gradient_field(cell, j, qp, dj) * gradient_field(cell, qp, dj);
                        }

                        assert(val == val);

                        auto dX = measure(cell, qp);
                        assert(dX == dX);

                        val += 1;
                        val *= fun(i, qp) * dX;
                        integral += val;
                    }

                    assert(integral == integral);
                    return integral;
                }

                Field gradient_field;
                Fun fun;
                Measure measure;
                const int n_qp;
                const int dim;
            };

            class Params : public Configurable {
            public:
                using Scalar = typename Traits<FirstLameParameter>::Scalar;
                void read(Input &) override {}
                Params() {}
            };

            LevelSet(const std::shared_ptr<FE> &fe, Params op = Params()) : Super(fe), op_(std::move(op)) {}

            inline int n_vars() const override { return 1; }

            inline bool is_matrix() const override { return true; }
            inline bool is_vector() const override { return true; }
            inline bool is_scalar() const override { return false; }
            inline bool is_linear() const override { return false; }

            inline std::string name() const override { return "LevelSet"; }

            virtual bool update(const std::shared_ptr<Field<FE>> &level_set_function) override {
                if (!Super::update(level_set_function)) {
                    return false;
                }

                assert(level_set_function);
                assert(level_set_function->is_coefficient());

                if (!level_set_function->is_coefficient()) {
                    Utopia::Abort("LevelSet::update, level_set_function must me in coefficient form!");
                }

                if (!gradient_field_) {
                    // Initialize gradient
                    gradient_field_ = std::make_shared<Gradient<FE>>(this->fe_ptr());
                }

                gradient_field_->init(*level_set_function);
                return true;
            }

            inline Op make_op() const {
                assert(level_set_function);
                auto grad_phi = gradient_field_->data();
                return Op(grad_phi, this->fe().grad(), this->fe().fun(), this->fe().measure());
            }

            inline GradientOp make_gradient_op() const {
                assert(level_set_function);
                auto grad_phi = gradient_field_->data();
                return GradientOp(grad_phi, this->fe().fun(), this->fe().measure());
            }

            bool apply(const DynRankView &x, DynRankView &y) override {
                UTOPIA_TRACE_REGION_BEGIN("LevelSet::apply");

                this->apply_vector_operator("LevelSet::apply", x, y, make_op());

                UTOPIA_TRACE_REGION_END("LevelSet::apply");
                return true;
            }

            bool assemble_matrix() override {
                UTOPIA_TRACE_REGION_BEGIN("LevelSet::assemble_matrix");

                this->loop_cell_test_trial("Assemble<LaplaceOperator>::assemble",
                                           op_and_store_cell_ij(this->matrix_data(), make_op()));

                UTOPIA_TRACE_REGION_END("LevelSet::assemble_matrix");
                return true;

                UTOPIA_TRACE_REGION_END("LevelSet::assemble_matrix");
                return true;
            }

            bool assemble_vector() override {
                UTOPIA_TRACE_REGION_BEGIN("LevelSet::assemble_vector");

                this->ensure_vector_accumulator();

                this->loop_cell_test("LevelSet::assemble_vector",
                                     op_and_store_cell_i(this->vector_data(), make_gradient_op()));

                UTOPIA_TRACE_REGION_END("LevelSet::assemble_vector");
                return true;
            }

            // NVCC_PRIVATE :
            Params op_;
            std::shared_ptr<Gradient<FE>> gradient_field_;
        };
    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_LEVELSET_HPP

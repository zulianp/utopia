#ifndef UTOPIA_KOKKOS_FEASSEMBLER_HPP
#define UTOPIA_KOKKOS_FEASSEMBLER_HPP

#include "utopia_Instance.hpp"
#include "utopia_Traits.hpp"

#include "utopia_FEAssembler.hpp"

#include "utopia_SimulationTime.hpp"
#include "utopia_fe_Core.hpp"

#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_Field.hpp"

#include "utopia_kokkos_TensorAccumulator.hpp"

#include "utopia_kokkos_SubdomainValue.hpp"

#include "utopia_kokkos_DiscretizationManager.hpp"

#include <memory>

namespace utopia {

    namespace kokkos {

        template <class View, class Op>
        class OpAndStoreCellIJ {
        public:
            template <typename... Args>
            UTOPIA_INLINE_FUNCTION OpAndStoreCellIJ(const View &data, Args &&...args)
                : op(std::forward<Args>(args)...), data(data) {}

            UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int i, const int j) const {
                data(cell, i, j) += op(cell, i, j);
            }

            Op op;
            View data;
        };

        template <class View, class Op>
        class BlockOpAndStoreCellIJ {
        public:
            template <typename... Args>
            UTOPIA_INLINE_FUNCTION BlockOpAndStoreCellIJ(const View &data, Args &&...args)
                : op(std::forward<Args>(args)...), data(data) {}

            UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int i, const int j) const {
                const int offset_i = i * op.dim();
                const int offset_j = j * op.dim();

                for (int sub_i = 0; sub_i < op.dim(); ++sub_i) {
                    const int dof_i = offset_i + sub_i;
                    for (int sub_j = 0; sub_j < op.dim(); ++sub_j) {
                        const int dof_j = offset_j + sub_j;
                        data(cell, dof_i, dof_j) += op(cell, i, j, sub_i, sub_j);
                    }
                }
            }

            Op op;
            View data;
        };

        template <class Op, class View, typename... Args>
        UTOPIA_INLINE_FUNCTION OpAndStoreCellIJ<View, Op> build_op_and_store_cell_ij(const View &data, Args &&...args) {
            return OpAndStoreCellIJ<View, Op>(data, std::forward<Args>(args)...);
        }

        template <class Op, class View>
        UTOPIA_INLINE_FUNCTION OpAndStoreCellIJ<View, Op> op_and_store_cell_ij(const View &data, Op op) {
            return OpAndStoreCellIJ<View, Op>(data, op);
        }

        template <class Op, class View>
        UTOPIA_INLINE_FUNCTION BlockOpAndStoreCellIJ<View, Op> block_op_and_store_cell_ij(const View &data, Op op) {
            return BlockOpAndStoreCellIJ<View, Op>(data, op);
        }

        template <class FunctionSpace,
                  typename FE,
                  typename MatrixView = DefaultView<typename FE::Scalar>,
                  typename VectorView = DefaultView<typename FE::Scalar>,
                  typename ScalarView = DefaultView<typename FE::Scalar>>
        class FEAssemblerWithManager : public FEAssembler<FunctionSpace> {
        public:
            using DiscretizationManager =
                utopia::kokkos::DiscretizationManager<FunctionSpace, FE, MatrixView, VectorView, ScalarView>;

            using Matrix = typename Traits<FunctionSpace>::Matrix;
            using Vector = typename Traits<FunctionSpace>::Vector;
            using Environment = utopia::Environment<FunctionSpace>;
            using Scalar = typename Traits<Vector>::Scalar;
            using SimulationTime = utopia::SimulationTime<Scalar>;

            virtual ~FEAssemblerWithManager() = default;

            void set_discretization_manager(const std::shared_ptr<DiscretizationManager> &dm) { dm_ = dm; }

            bool assemble(const Vector &x, Matrix &hessian, Vector &gradient) override {
                Utopia::Abort();
                return false;
            }
            bool assemble(const Vector &x, Matrix &jacobian) override {
                Utopia::Abort();
                return false;
            }
            bool assemble(const Vector &x, Vector &fun) override {
                Utopia::Abort();
                return false;
            }
            bool apply(const Vector &x, Vector &hessian_times_x) override {
                Utopia::Abort();
                return false;
            }

            // For linear only
            bool assemble(Matrix &jacobian) override {
                Utopia::Abort();
                return false;
            }
            bool assemble(Vector &fun) override {
                Utopia::Abort();
                return false;
            }

            void set_time(const std::shared_ptr<SimulationTime> &time) override { time_ = time; }

            void clear() override { Utopia::Abort(); }

            void set_environment(const std::shared_ptr<Environment> &env) override { Utopia::Abort(); }

            std::shared_ptr<Environment> environment() const override {
                Utopia::Abort();
                return nullptr;
            }

            void set_space(const std::shared_ptr<FunctionSpace> &space) override { Utopia::Abort(); }

            std::shared_ptr<FunctionSpace> space() const override {
                Utopia::Abort();
                return nullptr;
            }

            std::shared_ptr<SimulationTime> time() const { return time_; }

        private:
            std::shared_ptr<DiscretizationManager> dm_;
            std::shared_ptr<SimulationTime> time_;
        };

        // Handle case where no FunctionSpace is used i.e., FunctionSpace=void
        template <typename FE_, typename MatrixView_, typename VectorView_, typename ScalarView_>
        class FEAssemblerWithManager<void, FE_, MatrixView_, VectorView_, ScalarView_> : public FEAssembler<void> {
        public:
            virtual ~FEAssemblerWithManager() = default;
            bool assemble() { return false; }
        };

        template <class FunctionSpace_,
                  typename FE_,
                  typename MatrixView_ = DefaultView<typename FE_::Scalar>,
                  typename VectorView_ = DefaultView<typename FE_::Scalar>,
                  typename ScalarView_ = DefaultView<typename FE_::Scalar>>
        class FEAssembler : public Describable,
                            public FEAssemblerWithManager<FunctionSpace_, FE_, MatrixView_, VectorView_, ScalarView_> {
        public:
            using Super =
                utopia::kokkos::FEAssemblerWithManager<FunctionSpace_, FE_, MatrixView_, VectorView_, ScalarView_>;

            using FunctionSpace = FunctionSpace_;
            using FE = FE_;
            using MatrixView = MatrixView_;
            using VectorView = VectorView_;
            using ScalarView = ScalarView_;

            using Scalar = typename FE::Scalar;
            using SizeType = typename FE::SizeType;
            using ExecutionSpace = typename FE::ExecutionSpace;
            using MatrixAccumulator = utopia::kokkos::TensorAccumulator<MatrixView>;
            using VectorAccumulator = utopia::kokkos::TensorAccumulator<VectorView>;
            using ScalarAccumulator = utopia::kokkos::TensorAccumulator<ScalarView>;

            using CellTestTrialRange = typename FE::CellTestTrialRange;
            using CellTestRange = typename FE::CellTestRange;
            using CellRange = typename FE::CellRange;

            ////////////////////////////////////////////////////////////////////////////////////
            /// utopia::FEAssembler Types
            ////////////////////////////////////////////////////////////////////////////////////

            using Environment = typename Traits<FunctionSpace>::Environment;

            virtual ~FEAssembler() = default;

            virtual bool update(const std::shared_ptr<Field<FE>> &current_solution) {
                current_solution_ = current_solution;
                return true;
            }

            virtual bool apply(const VectorView &x, VectorView &y) {
                assert(false);
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                return false;
            }

            virtual int n_vars() const = 0;
            virtual std::string name() const = 0;

            void read(Input &) override {}

            FEAssembler(const std::shared_ptr<FE> &fe) : fe_(fe) { assert(fe); }

            inline void set_fe(const std::shared_ptr<FE> &fe) { fe_ = fe; }

            inline CellTestTrialRange cell_test_trial_range() const { return fe_->cell_test_trial_range(); }
            inline CellTestRange cell_test_range() const { return fe_->cell_test_range(); }
            inline CellRange cell_range() const { return fe_->cell_range(); }

            template <class CellFun>
            void loop_cell(const std::string &name, CellFun fun) const {
                Kokkos::parallel_for(name, cell_range(), fun);
            }

            template <class CellTestTrialFun>
            void loop_cell_test_trial(const std::string &name, CellTestTrialFun fun) const {
                Kokkos::parallel_for(name, cell_test_trial_range(), fun);
            }

            template <class CellTestFun>
            void loop_cell_test(const std::string &name, CellTestFun fun) const {
                Kokkos::parallel_for(name, cell_test_range(), fun);
            }

            template <class CellTestTrialOp>
            void apply_operator(const std::string &name, const VectorView &x, VectorView &y, CellTestTrialOp op) {
                assert(n_vars() == 1 && "Use apply_vector_operator");

                auto &fe = this->fe();

                const int n_shape_functions = fe.n_shape_functions();

                this->loop_cell_test(
                    name, UTOPIA_LAMBDA(const int &cell, const int &i) {
                        Scalar val = 0.0;
                        for (int j = 0; j < n_shape_functions; ++j) {
                            val += op(cell, i, j) * x(cell, j);
                        }

                        y(cell, i) += val;
                    });
            }

            template <class CellTestTrialSubIJOp>
            void apply_vector_operator(const std::string &name,
                                       const VectorView &x,
                                       VectorView &y,
                                       CellTestTrialSubIJOp op) {
                auto &fe = this->fe();

                const int n_shape_functions = fe.n_shape_functions();
                const int dim = op.dim();

                this->loop_cell_test(
                    name, UTOPIA_LAMBDA(const int &cell, const int &i) {
                        for (int j = 0; j < n_shape_functions; ++j) {
                            for (int sub_i = 0; sub_i < dim; ++sub_i) {
                                Scalar val = 0.0;
                                for (int sub_j = 0; sub_j < dim; ++sub_j) {
                                    val += op(cell, i, j, sub_i, sub_j) * x(cell, j * dim + sub_j);
                                }

                                y(cell, i * dim + sub_i) += val;
                            }
                        }
                    });
            }

            template <class CellTestOp>
            void apply_diagonal_operator(const std::string &name, const VectorView &x, VectorView &y, CellTestOp op) {
                assert(n_vars() == 1 && "IMPLEMENT ME");

                this->loop_cell_test(
                    name, UTOPIA_LAMBDA(const int &cell, const int &i) { y(cell, i) += op(cell, i) * x(cell, i); });
            }

            virtual bool is_matrix() const = 0;
            virtual bool is_vector() const = 0;
            virtual bool is_scalar() const = 0;

            virtual bool assemble_matrix() {
                assert(!is_matrix() && "IMPLEMENT ME IN SUB CLASS");
                return false;
            }

            virtual bool assemble_vector() {
                if (is_linear() && is_operator()) {
                    if (!this->current_solution()) {
                        assert(false);
                        return false;
                    }

                    this->ensure_vector_accumulator();

                    auto solution_field = this->current_solution();
                    assert(solution_field->is_coefficient());

                    auto data = this->vector_data();
                    return apply(solution_field->data(), data);
                } else {
                    assert(!is_vector() && "IMPLEMENT ME IN SUB CLASS");
                }

                return false;
            }
            virtual bool assemble_scalar() {
                assert(!is_scalar() && "IMPLEMENT ME IN SUB CLASS");
                return false;
            }

            virtual bool assemble() {
                if (is_scalar()) {
                    if (!assemble_scalar()) return false;
                }

                if (is_vector()) {
                    if (!assemble_vector()) return false;
                }

                if (is_matrix()) {
                    if (!assemble_matrix()) return false;
                }

                return true;
            }

            void scale(const kokkos::SubdomainValue<FE> &fun) {
                auto element_tags = fe().element_tags;

                if (is_matrix()) {
                    auto data = matrix_data();
                    scale_matrix(fun, data);
                }

                if (is_vector()) {
                    auto data = vector_data();
                    scale_vector(fun, data);
                }

                if (is_scalar()) {
                    auto data = scalar_data();
                    scale_scalar(fun, data);
                }
            }

            void scale_matrix(const kokkos::SubdomainValue<FE> &fun, MatrixView &matrix) const {
                auto element_tags = fe().element_tags();
                const int n_shape_functions = fe().n_shape_functions();

                loop_cell(
                    "FEAssembler::scale(mat)", UTOPIA_LAMBDA(int cell) {
                        auto tag = element_tags(cell);
                        auto val = fun.value(tag);

                        for (int i = 0; i < n_shape_functions; ++i) {
                            for (int j = 0; j < n_shape_functions; ++j) {
                                matrix(cell, i, j) *= val;
                            }
                        }
                    });
            }

            void scale_vector(const kokkos::SubdomainValue<FE> &fun, VectorView &vector) const {
                auto element_tags = fe().element_tags();
                const int n_shape_functions = fe().n_shape_functions();

                loop_cell(
                    "FEAssembler::scale_vector", UTOPIA_LAMBDA(int cell) {
                        auto val = fun.value(element_tags(cell));

                        for (int i = 0; i < n_shape_functions; ++i) {
                            vector(cell, i) *= val;
                        }
                    });
            }

            void scale_scalar(const kokkos::SubdomainValue<FE> &fun, ScalarView &scalar) const {
                auto element_tags = fe().element_tags();

                loop_cell(
                    "FEAssembler::scale_scalar", UTOPIA_LAMBDA(int cell) {
                        auto val = fun.value(element_tags(cell));
                        scalar(cell) *= val;
                    });
            }

            void ensure_accumulator() {
                if (is_scalar()) ensure_scalar_accumulator();
                if (is_vector()) ensure_vector_accumulator();
                if (is_matrix()) ensure_matrix_accumulator();
            }

            inline std::shared_ptr<MatrixAccumulator> matrix_accumulator() const { return matrix_accumulator_; }
            inline std::shared_ptr<VectorAccumulator> vector_accumulator() const { return vector_accumulator_; }
            inline std::shared_ptr<ScalarAccumulator> scalar_accumulator() const { return scalar_accumulator_; }

            virtual void set_matrix_accumulator(const std::shared_ptr<MatrixAccumulator> &matrix_accumulator) {
                matrix_accumulator_ = matrix_accumulator;
            }

            virtual void set_vector_accumulator(const std::shared_ptr<VectorAccumulator> &vector_accumulator) {
                vector_accumulator_ = vector_accumulator;
            }

            virtual void set_scalar_accumulator(const std::shared_ptr<ScalarAccumulator> &scalar_accumulator) {
                scalar_accumulator_ = scalar_accumulator;
            }

            std::shared_ptr<FE> fe_ptr() { return fe_; }
            FE &fe() { return *fe_; }
            const FE &fe() const { return *fe_; }

            inline MatrixView matrix_data() {
                assert(matrix_accumulator_);
                if (!matrix_accumulator_) {
                    return MatrixView();
                } else {
                    return matrix_accumulator_->data();
                }
            }

            inline VectorView vector_data() {
                assert(vector_accumulator_);
                if (!vector_accumulator_) {
                    return VectorView();
                } else {
                    return vector_accumulator_->data();
                }
            }

            inline ScalarView scalar_data() {
                assert(scalar_accumulator_);
                if (!scalar_accumulator_) {
                    return ScalarView();
                } else {
                    return scalar_accumulator_->data();
                }
            }

            void describe(std::ostream &os) const override {
                const SizeType n_cells = fe_->n_cells();
                const int n_shape_functions = fe_->n_shape_functions();
                os << "n_cells: " << n_cells << ", n_shape_functions: " << n_shape_functions << "\n";

                if (matrix_accumulator_) {
                    matrix_accumulator_->describe(os);
                }

                if (vector_accumulator_) {
                    vector_accumulator_->describe(os);
                }

                if (scalar_accumulator_) {
                    scalar_accumulator_->describe(os);
                }
            }

            virtual void ensure_matrix_accumulator() {
                if (!matrix_accumulator_) {
                    matrix_accumulator_ = std::make_shared<MatrixAccumulator>();
                    matrix_accumulator_->init_matrix(*fe_, n_vars());
                } else {
                    matrix_accumulator_->prepare();
                }
            }

            virtual void ensure_vector_accumulator() {
                if (!vector_accumulator_) {
                    vector_accumulator_ = std::make_shared<VectorAccumulator>();
                    vector_accumulator_->init_vector(*fe_, n_vars());
                } else {
                    vector_accumulator_->prepare();
                }
            }

            virtual void ensure_scalar_accumulator() {
                if (!scalar_accumulator_) {
                    scalar_accumulator_ = std::make_shared<ScalarAccumulator>();
                    scalar_accumulator_->init_scalar(*fe_, n_vars());
                } else {
                    scalar_accumulator_->prepare();
                }
            }

            inline std::shared_ptr<Field<FE>> current_solution() { return current_solution_; }

            //////////////////////////////////////////////////////////////////////////////

            bool is_operator() const override { return false; }
            bool is_linear() const override { return true; }

        private:
            std::shared_ptr<FE> fe_;
            std::shared_ptr<MatrixAccumulator> matrix_accumulator_;
            std::shared_ptr<VectorAccumulator> vector_accumulator_;
            std::shared_ptr<ScalarAccumulator> scalar_accumulator_;
            std::shared_ptr<Field<FE>> current_solution_;
        };

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_FEASSEMBLER_HPP

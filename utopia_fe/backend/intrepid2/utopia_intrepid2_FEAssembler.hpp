#ifndef UTOPIA_INTREPID2_FEASSEMBLER_HPP
#define UTOPIA_INTREPID2_FEASSEMBLER_HPP

#include "utopia_fe_Core.hpp"
#include "utopia_intrepid2_Commons.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_ForwardDeclarations.hpp"

#include <memory>

namespace utopia {
    namespace intrepid2 {
        template <class Operator, typename Scalar = UScalar>
        class Assemble {};

        template <class View, class Op>
        class OpAndStoreCellIJ {
        public:
            template <typename... Args>
            UTOPIA_INLINE_FUNCTION OpAndStoreCellIJ(const View &data, Args &&... args)
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
            UTOPIA_INLINE_FUNCTION BlockOpAndStoreCellIJ(const View &data, Args &&... args)
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
        UTOPIA_INLINE_FUNCTION OpAndStoreCellIJ<View, Op> build_op_and_store_cell_ij(const View &data,
                                                                                     Args &&... args) {
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

        template <typename Scalar>
        class FEAssembler : public Describable, public Configurable {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using ExecutionSpace = typename FE::ExecutionSpace;

            using CellTestTrialRange = typename FE::CellTestTrialRange;
            using CellTestRange = typename FE::CellTestRange;
            using CellRange = typename FE::CellRange;

            virtual ~FEAssembler() = default;

            virtual bool update(const std::shared_ptr<Field<Scalar>> &current_solution) {
                current_solution_ = current_solution;
                return true;
            }

            virtual bool apply(const DynRankView &x, DynRankView &y) {
                assert(false);
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                return false;
            }

            virtual int n_vars() const = 0;
            virtual std::string name() const = 0;

            void read(Input &) override {}

            FEAssembler(const std::shared_ptr<FE> &fe) : fe_(fe) { assert(fe); }

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
            void apply_operator(const std::string &name, const DynRankView &x, DynRankView &y, CellTestTrialOp op) {
                assert(n_vars() == 1 && "IMPLEMENT ME");

                auto &fe = this->fe();

                const int num_fields = fe.num_fields();

                this->loop_cell_test(
                    name, UTOPIA_LAMBDA(const int &cell, const int &i) {
                        Scalar val = 0.0;
                        for (int j = 0; j < num_fields; ++j) {
                            val += op(cell, i, j) * x(cell, j);
                        }

                        y(cell, i) += val;
                    });
            }

            template <class CellTestTrialSubIJOp>
            void apply_vector_operator(const std::string &name,
                                       const DynRankView &x,
                                       DynRankView &y,
                                       CellTestTrialSubIJOp op) {
                assert(n_vars() == 1 && "IMPLEMENT ME");

                auto &fe = this->fe();

                const int num_fields = fe.num_fields();
                const int dim = op.dim();

                this->loop_cell_test(
                    name, UTOPIA_LAMBDA(const int &cell, const int &i) {
                        for (int j = 0; j < num_fields; ++j) {
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
            void apply_diagonal_operator(const std::string &name, const DynRankView &x, DynRankView &y, CellTestOp op) {
                assert(n_vars() == 1 && "IMPLEMENT ME");

                this->loop_cell_test(
                    name, UTOPIA_LAMBDA(const int &cell, const int &i) { y(cell, i) += op(cell, i) * x(cell, i); });
            }

            virtual bool is_matrix() const = 0;
            virtual bool is_vector() const = 0;
            virtual bool is_scalar() const = 0;
            virtual bool is_operator() const { return false; }

            virtual bool is_linear() const { return true; }

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

            void scale(const intrepid2::SubdomainValue<Scalar> &fun) {
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

            void scale_matrix(const intrepid2::SubdomainValue<Scalar> &fun, DynRankView &matrix) const {
                auto element_tags = fe().element_tags;
                const int num_fields = fe().num_fields();

                loop_cell(
                    "FEAssembler::scale(mat)", UTOPIA_LAMBDA(int cell) {
                        auto val = fun.value(element_tags(cell));

                        for (int i = 0; i < num_fields; ++i) {
                            for (int j = 0; j < num_fields; ++j) {
                                matrix(cell, i, j) *= val;
                            }
                        }
                    });
            }

            void scale_vector(const intrepid2::SubdomainValue<Scalar> &fun, DynRankView &vector) const {
                auto element_tags = fe().element_tags;
                const int num_fields = fe().num_fields();

                loop_cell(
                    "FEAssembler::scale_vector", UTOPIA_LAMBDA(int cell) {
                        auto val = fun.value(element_tags(cell));

                        for (int i = 0; i < num_fields; ++i) {
                            vector(cell, i) *= val;
                        }
                    });
            }

            void scale_scalar(const intrepid2::SubdomainValue<Scalar> &fun, DynRankView &scalar) const {
                auto element_tags = fe().element_tags;

                loop_cell(
                    "FEAssembler::scale_scalar", UTOPIA_LAMBDA(int cell) {
                        auto val = fun.value(element_tags(cell));
                        scalar(cell) *= val;
                    });
            }

            class TensorAccumulator : public Describable {
            public:
                DynRankView &data() { return data_; }

                void init_scalar(const FE &fe, int n_vars) { data_ = DynRankView("scalars", fe.num_cells(), n_vars); }

                void init_matrix(const FE &fe, int n_vars) {
                    const int num_fields = fe.num_fields();
                    const int n_dofs = num_fields * n_vars;
                    data_ = DynRankView("matrices", fe.num_cells(), n_dofs, n_dofs);
                }

                void init_vector(FE &fe, int n_vars) {
                    const int num_fields = fe.num_fields();
                    const int n_dofs = num_fields * n_vars;
                    data_ = DynRankView("vectors", fe.num_cells(), n_dofs);
                }

                inline AssemblyMode mode() const { return mode_; }
                inline void set_mode(AssemblyMode mode) const { mode_ = mode; }

                bool is_compatible(const FE &fe) const { return fe.num_cells() <= data_.extent(0); }

                void prepare() {
                    if (mode_ == OVERWRITE_MODE) {
                        zero();
                    }
                }

                void zero() { utopia::intrepid2::fill(data_, 0.0); }

                void describe(std::ostream &os) const override {
                    const SizeType num_cells = data_.extent(0);

                    if (data_.rank() == 3) {
                        const int n_dofs_i = data_.extent(1);

                        for (SizeType c = 0; c < num_cells; ++c) {
                            os << c << ")\n";
                            for (SizeType i = 0; i < n_dofs_i; ++i) {
                                os << data_(c, i) << " ";

                                os << '\n';
                            }

                            os << '\n';
                        }

                    } else if (data_.rank() == 4) {
                        const int n_dofs_i = data_.extent(1);
                        const int n_dofs_j = data_.extent(2);

                        for (SizeType c = 0; c < num_cells; ++c) {
                            os << c << ")\n";
                            for (SizeType i = 0; i < n_dofs_i; ++i) {
                                for (SizeType j = 0; j < n_dofs_j; ++j) {
                                    os << data_(c, i, j) << " ";
                                }

                                os << '\n';
                            }

                            os << '\n';
                        }
                    }
                }

            private:
                DynRankView data_;
                AssemblyMode mode_{ADD_MODE};
            };

            void ensure_accumulator() {
                if (is_scalar()) ensure_scalar_accumulator();
                if (is_vector()) ensure_vector_accumulator();
                if (is_matrix()) ensure_matrix_accumulator();
            }

            inline std::shared_ptr<TensorAccumulator> matrix_accumulator() { return matrix_accumulator_; }
            inline std::shared_ptr<TensorAccumulator> matrix_accumulator() const { return matrix_accumulator_; }

            inline std::shared_ptr<TensorAccumulator> vector_accumulator() { return vector_accumulator_; }
            inline std::shared_ptr<TensorAccumulator> vector_accumulator() const { return vector_accumulator_; }

            inline std::shared_ptr<TensorAccumulator> scalar_accumulator() { return scalar_accumulator_; }
            inline std::shared_ptr<TensorAccumulator> scalar_accumulator() const { return scalar_accumulator_; }

            void set_matrix_accumulator(const std::shared_ptr<TensorAccumulator> &matrix_accumulator) {
                matrix_accumulator_ = matrix_accumulator;
            }

            void set_vector_accumulator(const std::shared_ptr<TensorAccumulator> &vector_accumulator) {
                vector_accumulator_ = vector_accumulator;
            }

            void set_scalar_accumulator(const std::shared_ptr<TensorAccumulator> &scalar_accumulator) {
                scalar_accumulator_ = scalar_accumulator;
            }

            std::shared_ptr<FE> fe_ptr() { return fe_; }
            FE &fe() { return *fe_; }
            const FE &fe() const { return *fe_; }

            inline DynRankView matrix_data() {
                assert(matrix_accumulator_);
                if (!matrix_accumulator_) {
                    return DynRankView();
                } else {
                    return matrix_accumulator_->data();
                }
            }

            inline DynRankView vector_data() {
                assert(vector_accumulator_);
                if (!vector_accumulator_) {
                    return DynRankView();
                } else {
                    return vector_accumulator_->data();
                }
            }

            inline DynRankView scalar_data() {
                assert(scalar_accumulator_);
                if (!scalar_accumulator_) {
                    return DynRankView();
                } else {
                    return scalar_accumulator_->data();
                }
            }

            void describe(std::ostream &os) const override {
                const SizeType num_cells = fe_->num_cells();
                const int num_fields = fe_->num_fields();
                os << "num_cells: " << num_cells << ", num_fields: " << num_fields << "\n";

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

            void ensure_matrix_accumulator() {
                if (!matrix_accumulator_) {
                    matrix_accumulator_ = std::make_shared<TensorAccumulator>();
                    matrix_accumulator_->init_matrix(*fe_, n_vars());
                } else {
                    matrix_accumulator_->prepare();
                }
            }

            void ensure_vector_accumulator() {
                if (!vector_accumulator_) {
                    vector_accumulator_ = std::make_shared<TensorAccumulator>();
                    vector_accumulator_->init_vector(*fe_, n_vars());
                } else {
                    vector_accumulator_->prepare();
                }
            }

            void ensure_scalar_accumulator() {
                if (!scalar_accumulator_) {
                    scalar_accumulator_ = std::make_shared<TensorAccumulator>();
                    scalar_accumulator_->init_scalar(*fe_, n_vars());
                } else {
                    scalar_accumulator_->prepare();
                }
            }

            inline std::shared_ptr<Field<Scalar>> current_solution() { return current_solution_; }

        private:
            std::shared_ptr<FE> fe_;
            std::shared_ptr<TensorAccumulator> matrix_accumulator_;
            std::shared_ptr<TensorAccumulator> vector_accumulator_;
            std::shared_ptr<TensorAccumulator> scalar_accumulator_;
            std::shared_ptr<Field<Scalar>> current_solution_;
        };

    }  // namespace intrepid2

    template <typename Scalar, class Op>
    class AssembleTraits<intrepid2::FE<Scalar>, Op> {
    public:
        using Type = utopia::intrepid2::Assemble<Op, Scalar>;
    };
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FEASSEMBLER_HPP
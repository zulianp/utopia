#ifndef UTOPIA_INTREPID2_FEASSEMBLER_HPP
#define UTOPIA_INTREPID2_FEASSEMBLER_HPP

#include "utopia_fe_Core.hpp"
#include "utopia_intrepid2_Commons.hpp"

#include "utopia_intrepid2_FE.hpp"
#include "utopia_intrepid2_ForwardDeclarations.hpp"

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

            UTOPIA_INLINE_FUNCTION void operator()(const int &cell, const int &i, const int &j) const {
                data(cell, i, j) += op(cell, i, j);
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

        template <typename Scalar>
        class FEAssembler : public Describable, public Configurable {
        public:
            using FE = utopia::intrepid2::FE<Scalar>;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using FunctionSpaceTools = typename FE::FunctionSpaceTools;
            using ExecutionSpace = typename FE::ExecutionSpace;

            virtual ~FEAssembler() = default;
            virtual bool assemble() = 0;

            virtual bool apply(const DynRankView &x, DynRankView &y) {
                assert(false);
                UTOPIA_UNUSED(x);
                UTOPIA_UNUSED(y);
                return false;
            }
            virtual int n_vars() const = 0;
            virtual std::string name() const = 0;
            virtual int rank() const = 0;

            void read(Input &) override {}

            FEAssembler(const std::shared_ptr<FE> &fe) : fe_(fe) { assert(fe); }

            Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace> cell_test_trial_range() const {
                int num_cells = fe_->num_cells();
                int num_fields = fe_->num_fields();

                return Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0},
                                                                              {num_cells, num_fields, num_fields});
            }

            Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace> cell_test_range() const {
                int num_cells = fe_->num_cells();
                int num_fields = fe_->num_fields();

                return Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>({0, 0}, {num_cells, num_fields});
            }

            Kokkos::RangePolicy<ExecutionSpace> cell_range() const {
                int num_cells = fe_->num_cells();

                return Kokkos::RangePolicy<ExecutionSpace>(0, num_cells);
            }

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

            template <class CellTestOp>
            void apply_diagonal_operator(const std::string &name, const DynRankView &x, DynRankView &y, CellTestOp op) {
                assert(n_vars() == 1 && "IMPLEMENT ME");

                this->loop_cell_test(
                    name, UTOPIA_LAMBDA(const int &cell, const int &i) { y(cell, i) += op(cell, i) * x(cell, i); });
            }

            bool is_matrix() const {
                if (!accumulator()) return false;
                return accumulator()->data().rank() == 4;
            }

            bool is_vector() const {
                if (!accumulator()) return false;
                return accumulator()->data().rank() == 3;
            }

            bool is_scalar() const {
                if (!accumulator()) return false;
                return accumulator()->data().rank() == 2;
            }

            void scale(const intrepid2::SubdomainValue<Scalar> &fun) {
                auto element_tags = fe().element_tags;
                const int num_fields = fe().num_fields();
                auto data = accumulator()->data();

                if (is_matrix()) {
                    loop_cell(
                        "FEAssembler::scale(mat)", UTOPIA_LAMBDA(int cell) {
                            auto val = fun.value(element_tags(cell));

                            for (int i = 0; i < num_fields; ++i) {
                                for (int j = 0; j < num_fields; ++j) {
                                    data(cell, i, j) *= val;
                                }
                            }
                        });
                } else if (is_vector()) {
                    scale_vector(fun, data);
                } else if (is_scalar()) {
                    loop_cell(
                        "FEAssembler::scale(scalar)", UTOPIA_LAMBDA(int cell) {
                            auto val = fun.value(element_tags(cell));

                            for (int i = 0; i < num_fields; ++i) {
                                data(cell) *= val;
                            }
                        });
                } else {
                    assert(false);
                    Utopia::Abort("Called FEAssembler::scale on invalid accumulator!");
                }
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

            class TensorAccumulator {
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

            private:
                DynRankView data_;
                AssemblyMode mode_{ADD_MODE};
            };

            void ensure_accumulator() {
                switch (rank()) {
                    case 0: {
                        ensure_scalar_accumulator();
                        break;
                    }
                    case 1: {
                        ensure_vec_accumulator();
                        break;
                    }
                    case 2: {
                        ensure_mat_accumulator();
                        break;
                    }
                    default: {
                        assert(false && "INVALID RANK");
                        Utopia::Abort();
                    }
                }
            }

            inline std::shared_ptr<TensorAccumulator> accumulator() { return accumulator_; }
            inline std::shared_ptr<TensorAccumulator> accumulator() const { return accumulator_; }

            void set_accumulator(const std::shared_ptr<TensorAccumulator> &accumulator) { accumulator_ = accumulator; }

            std::shared_ptr<FE> fe_ptr() { return fe_; }
            FE &fe() { return *fe_; }
            const FE &fe() const { return *fe_; }

            inline DynRankView data() {
                assert(accumulator_);
                if (!accumulator_) {
                    return DynRankView();
                } else {
                    return accumulator_->data();
                }
            }

            void describe(std::ostream &os) const override {
                if (accumulator_) {
                    auto &data = accumulator_->data();

                    const SizeType num_cells = fe_->num_cells();
                    const int num_fields = fe_->num_fields();
                    const int n_dofs_i = data.extent(1);

                    os << "num_cells: " << num_cells << ", num_fields: " << num_fields << "\n";

                    if (data.rank() == 3) {
                        for (SizeType c = 0; c < num_cells; ++c) {
                            os << c << ")\n";
                            for (SizeType i = 0; i < n_dofs_i; ++i) {
                                os << data(c, i) << " ";

                                os << '\n';
                            }

                            os << '\n';
                        }

                    } else if (data.rank() == 4) {
                        const int n_dofs_j = data.extent(2);

                        for (SizeType c = 0; c < num_cells; ++c) {
                            os << c << ")\n";
                            for (SizeType i = 0; i < n_dofs_i; ++i) {
                                for (SizeType j = 0; j < n_dofs_j; ++j) {
                                    os << data(c, i, j) << " ";
                                }

                                os << '\n';
                            }

                            os << '\n';
                        }
                    }
                }
            }

        private:
            std::shared_ptr<FE> fe_;
            std::shared_ptr<TensorAccumulator> accumulator_;

        protected:
            void ensure_mat_accumulator() {
                if (!accumulator_) {
                    accumulator_ = std::make_shared<TensorAccumulator>();
                    accumulator_->init_matrix(*fe_, n_vars());
                } else {
                    accumulator_->prepare();
                }
            }

            void ensure_vec_accumulator() {
                if (!accumulator_) {
                    accumulator_ = std::make_shared<TensorAccumulator>();
                    accumulator_->init_vector(*fe_, n_vars());
                } else {
                    accumulator_->prepare();
                }
            }

            void ensure_scalar_accumulator() {
                if (!accumulator_) {
                    accumulator_ = std::make_shared<TensorAccumulator>();
                    accumulator_->init_scalar(*fe_, n_vars());
                } else {
                    accumulator_->prepare();
                }
            }
        };

    }  // namespace intrepid2

    template <typename Scalar, class Op>
    class AssembleTraits<intrepid2::FE<Scalar>, Op> {
    public:
        using Type = utopia::intrepid2::Assemble<Op, Scalar>;
    };
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FEASSEMBLER_HPP
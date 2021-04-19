#ifndef UTOPIA_INTREPID2_FEASSEMBLER_HPP
#define UTOPIA_INTREPID2_FEASSEMBLER_HPP

#include "utopia_fe_Core.hpp"
#include "utopia_intrepid2_Commons.hpp"

#include "utopia_intrepid2_FE.hpp"

namespace utopia {
    namespace intrepid2 {
        template <class Operator, typename Scalar = UScalar>
        class Assemble {};

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
            virtual int n_vars() const = 0;
            virtual std::string name() const = 0;
            virtual int rank() const = 0;

            void read(Input &) override {}

            FEAssembler(const std::shared_ptr<FE> &fe) : fe_(fe) { assert(fe); }

            Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace> mat_range() {
                int num_cells = fe_->num_cells();
                int num_fields = fe_->num_fields();

                return Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>({0, 0, 0},
                                                                              {num_cells, num_fields, num_fields});
            }

            Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace> vec_range() {
                int num_cells = fe_->num_cells();
                int num_fields = fe_->num_fields();

                return Kokkos::MDRangePolicy<Kokkos::Rank<2>, ExecutionSpace>({0, 0}, {num_cells, num_fields});
            }

            template <class CellIJFun>
            void mat_integrate(const std::string &name, CellIJFun fun) {
                Kokkos::parallel_for(name, mat_range(), fun);
            }

            template <class CellIFun>
            void vec_integrate(const std::string &name, CellIFun fun) {
                Kokkos::parallel_for(name, vec_range(), fun);
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
}  // namespace utopia

#endif  // UTOPIA_INTREPID2_FEASSEMBLER_HPP
#ifndef UTOPIA_KOKKOS_TENSOR_FIELD_HPP
#define UTOPIA_KOKKOS_TENSOR_FIELD_HPP

#include "utopia_Traits.hpp"
#include "utopia_Views.hpp"

#include "utopia_DeviceOperations.hpp"

#include "utopia_kokkos_Commons.hpp"
#include "utopia_kokkos_FE.hpp"
#include "utopia_kokkos_FEAssembler.hpp"
#include "utopia_kokkos_Field.hpp"

namespace utopia {
    namespace kokkos {

        template <typename FE_>
        class QPTensorField : public QPField<FE_> {
        public:
            using FE = FE_;
            using Super = utopia::kokkos::QPField<FE>;
            using Scalar = typename FE::Scalar;
            using SizeType = typename FE::SizeType;
            using DynRankView = typename FE::DynRankView;
            using ExecutionSpace = typename FE::ExecutionSpace;

            using Rank1Range = ::Kokkos::MDRangePolicy<::Kokkos::Rank<3>, ExecutionSpace>;
            using Rank2Range = ::Kokkos::MDRangePolicy<::Kokkos::Rank<4>, ExecutionSpace>;

            QPTensorField(const std::shared_ptr<FE> &fe, const std::string &name = "QPTensorField") : Super(fe) {
                this->set_name(name);
            }

            void set_tensor_size(const int rows, const int cols) {
                rows_ = rows;
                cols_ = cols;
                Super::set_tensor_size(rows * cols);
            }

            void ensure_field() override {
                const SizeType e0 = this->data().extent(0);
                const SizeType e1 = this->data().extent(1);
                const SizeType e2 = this->data().extent(2);

                if (this->data().rank() != 3 || e0 < this->fe()->n_cells() || e1 < this->fe()->n_quad_points() ||
                    e2 < this->tensor_size()) {
                    this->data() = DynRankView(
                        this->name(), this->fe()->n_cells(), this->fe()->n_quad_points(), this->tensor_size());
                } else {
                    utopia::kokkos::fill(this->data(), 0.0);
                }
            }

            inline int rank() const { return (rows_ == 1 || cols_ == 1) ? 1 : 2; }
            bool is_coefficient() const override { return false; }

            inline int rows() const { return rows_; }
            inline int cols() const { return cols_; }

            class AddIdentity {
            public:
                UTOPIA_INLINE_FUNCTION AddIdentity(const int &rows, const int &cols, const DynRankView &field)
                    : rows(rows), cols(cols), field(field) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(const int cell, const int qp, const int r, const int c) const {
                    return field(cell, qp, r * cols + c) + (r == c);
                }

                const int rows, cols;
                const DynRankView field;
            };

            class AddIdentityAndStore {
            public:
                AddIdentityAndStore(const int &rows, const int &cols, DynRankView &field) : op_(rows, cols, field) {}

                UTOPIA_INLINE_FUNCTION void operator()(const int cell, const int qp) const {
                    for (int r = 0; r < op_.rows; ++r) {
                        for (int c = 0; c < op_.cols; ++c) {
                            op_.field(cell, qp, r * op_.cols + c) = op_(cell, qp, r, c);
                        }
                    }
                }

                AddIdentity op_;
            };

            void add_identity() {
                auto data = this->data();
                assert(rows() == cols());
                ::Kokkos::parallel_for(this->fe()->cell_qp_range(), AddIdentityAndStore(rows(), cols(), data));
            }

            void describe(std::ostream &) const override {
                const int rows = this->rows();
                const int cols = this->cols();

                auto data = this->data();
                ::Kokkos::parallel_for(
                    this->fe()->cell_qp_range(), UTOPIA_LAMBDA(int cell, int qp) {
                        if (qp == 0) {
                            printf("\n");
                        }

                        for (int i = 0; i < rows; ++i) {
                            for (int j = 0; j < cols; ++j) {
                                printf("%g ", data(cell, qp, i * cols + j));
                            }

                            printf("\n");
                        }

                        printf("\n");
                    });
            }

            class EigenvaluesOpAndStore {
            public:
                EigenvaluesOpAndStore(const int n, const DynRankView &data, DynRankView &result)
                    : n(n), data(data), result(result) {}

                UTOPIA_INLINE_FUNCTION void operator()(int cell, int qp) const {
                    switch (n) {
                        case 1: {
                            result(cell, qp, 0) = data(cell, qp, 0);
                            break;
                        }

                        case 2: {
                            StaticVector<Scalar, 2> e;
                            StaticMatrix<Scalar, 2, 2> mat;

                            mat(0, 0) = data(cell, qp, 0);
                            mat(0, 1) = data(cell, qp, 1);
                            mat(1, 0) = data(cell, qp, 2);
                            mat(1, 1) = data(cell, qp, 3);

                            utopia::eig(mat, e);

                            result(cell, qp, 0) = e[0];
                            result(cell, qp, 1) = e[1];
                            break;
                        }

                        case 3: {
                            StaticVector<Scalar, 3> e;
                            StaticMatrix<Scalar, 3, 3> mat;

                            mat(0, 0) = data(cell, qp, 0);
                            mat(0, 1) = data(cell, qp, 1);
                            mat(0, 2) = data(cell, qp, 2);

                            mat(1, 0) = data(cell, qp, 3);
                            mat(1, 1) = data(cell, qp, 4);
                            mat(1, 2) = data(cell, qp, 5);

                            mat(2, 0) = data(cell, qp, 6);
                            mat(2, 1) = data(cell, qp, 7);
                            mat(2, 2) = data(cell, qp, 8);

                            utopia::eig(mat, e);

                            result(cell, qp, 0) = e[0];
                            result(cell, qp, 1) = e[1];
                            result(cell, qp, 2) = e[2];
                            break;
                        }

                        default: {
                            assert(false);
                            break;
                        }
                    }
                }

                int n;
                DynRankView data;
                DynRankView result;
            };

            class DetOp {
            public:
                DetOp(const int n, const DynRankView &data) : n(n), data(data) {}

                UTOPIA_INLINE_FUNCTION Scalar operator()(int cell, int qp) const {
                    Scalar d = 0;

                    switch (n) {
                        case 1: {
                            d = data(cell, qp, 0);
                            break;
                        }

                        case 2: {
                            StaticMatrix<Scalar, 2, 2> mat;
                            mat(0, 0) = data(cell, qp, 0);
                            mat(0, 1) = data(cell, qp, 1);
                            mat(1, 0) = data(cell, qp, 2);
                            mat(1, 1) = data(cell, qp, 3);
                            d = utopia::det(mat);
                            break;
                        }

                        case 3: {
                            StaticMatrix<Scalar, 3, 3> mat;
                            mat(0, 0) = data(cell, qp, 0);
                            mat(0, 1) = data(cell, qp, 1);
                            mat(0, 2) = data(cell, qp, 2);

                            mat(1, 0) = data(cell, qp, 3);
                            mat(1, 1) = data(cell, qp, 4);
                            mat(1, 2) = data(cell, qp, 5);

                            mat(2, 0) = data(cell, qp, 6);
                            mat(2, 1) = data(cell, qp, 7);
                            mat(2, 2) = data(cell, qp, 8);

                            d = utopia::det(mat);
                            break;
                        }

                        default: {
                            assert(false);
                            break;
                        }
                    }

                    return d;
                }

                int n;
                DynRankView data;
            };

            class CountNonPositiveDets {
            public:
                CountNonPositiveDets(const int n, const DynRankView &matrices)
                    : op_(n, matrices), n_qp(matrices.extent(1)) {}

                UTOPIA_INLINE_FUNCTION int apply(int cell, int qp) const {
                    Scalar J = op_(cell, qp);
                    bool is_non_positive = J <= 0.0;
                    assert(!is_non_positive);
                    return is_non_positive;
                }

                KOKKOS_INLINE_FUNCTION void operator()(const int &cell, int &val) const {
                    for (int qp = 0; qp < n_qp; ++qp) {
                        val += apply(cell, qp);
                    }
                }

                UTOPIA_INLINE_FUNCTION void init(int &val) const { val = 0; }

#define UTOPIA_KOKKOS_JOIN_VOLATILE volatile

#ifdef KOKKOS_VERSION_GREATER_EQUAL
#if KOKKOS_VERSION_GREATER_EQUAL(4, 0, 1)
#define UTOPIA_KOKKOS_JOIN_VOLATILE
#endif
#endif

                KOKKOS_INLINE_FUNCTION void join(UTOPIA_KOKKOS_JOIN_VOLATILE int &val,
                                                 const UTOPIA_KOKKOS_JOIN_VOLATILE int &other) const

                {
                    // Kokkos forces us to have the input values being declared volatile. Hence we need to make copies
                    // for the reduction operations
                    const int tmp1 = val, tmp2 = other;
                    val = tmp1 + tmp2;
                }

                DetOp op_;
                int n_qp;
            };

            bool check_dets_are_positive() const {
                if (rows() != cols()) {
                    assert(false);
                    Utopia::Abort("Trying to call det on non-square QPTensorField");
                }

                CountNonPositiveDets op(this->rows(), this->data());

                int count = 0;
                ::Kokkos::parallel_reduce(this->fe()->cell_range(), op, count);
                assert(count == 0);
                return count == 0;
            }

            void det(DynRankView &result) const {
                if (rows() != cols()) {
                    assert(false);
                    Utopia::Abort("Trying to call det on non-square QPTensorField");
                }

                SizeType n_cells = this->fe()->n_cells();
                SizeType n_quad_points = this->fe()->n_quad_points();

                if (SizeType(result.extent(0)) != n_cells || SizeType(result.extent(1)) != n_quad_points ||
                    SizeType(result.extent(2)) != 1) {
                    result = DynRankView("det(QPTensorField)", n_cells, n_quad_points, 1);
                }

                const int n = this->rows();

                auto data = this->data();

                DetOp op(n, data);
                ::Kokkos::parallel_for(
                    this->fe()->cell_qp_range(),
                    UTOPIA_LAMBDA(int cell, int qp) { result(cell, qp, 0) = op(cell, qp); });
            }

            void eig(QPField<FE> &field) {
                field.set_tensor_size(rows());
                eig(field.data());
            }

            void eig(DynRankView &result) const {
                if (rows() != cols()) {
                    assert(false);
                    Utopia::Abort("Trying to call eig on non-square QPTensorField");
                }

                const SizeType n_cells = this->fe()->n_cells();
                const SizeType n_quad_points = this->fe()->n_quad_points();
                const SizeType n = this->rows();

                if (SizeType(result.extent(0)) != n_cells || SizeType(result.extent(1)) != n_quad_points ||
                    SizeType(result.extent(2)) != n) {
                    result = DynRankView("eig(QPTensorField)", n_cells, n_quad_points, n);
                }

                auto data = this->data();

                EigenvaluesOpAndStore op(n, data, result);
                ::Kokkos::parallel_for(
                    this->fe()->cell_qp_range(), UTOPIA_LAMBDA(int cell, int qp) { op(cell, qp); });
            }

            inline Rank1Range rank1_range() {
                int n_cells = this->fe()->n_cells();
                int n_quad_points = this->fe()->n_quad_points();
                return Kokkos::MDRangePolicy<Kokkos::Rank<3>, ExecutionSpace>(
                    {0, 0, 0}, {n_cells, n_quad_points, this->tensor_size()});
            }

            inline Rank2Range rank2_range() {
                int n_cells = this->fe()->n_cells();
                int n_quad_points = this->fe()->n_quad_points();
                return Kokkos::MDRangePolicy<Kokkos::Rank<4>, ExecutionSpace>({0, 0, 0, 0},
                                                                              {n_cells, n_quad_points, rows(), cols()});
            }

        private:
            int rows_{1}, cols_{1};
        };

        template <typename FE>
        inline QPField<FE> det(const QPTensorField<FE> &grad) {
            QPField<FE> ret(grad.fe());
            ret.set_name("Determinant");
            grad.det(ret.data());
            return ret;
        }

    }  // namespace kokkos
}  // namespace utopia

#endif  // UTOPIA_KOKKOS_TENSOR_FIELD_HPP
#include "utopia_moonolith_FETransfer.hpp"
#include "utopia_ExtractComponent.hpp"

#include "utopia_ElementWisePseudoInverse.hpp"
#include "utopia_Options.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_moonolith_ConvertTensor.hpp"

// All Moonolith includes
#include "moonolith_communicator.hpp"
#include "moonolith_par_l2_assembler_no_covering.hpp"
#include "moonolith_par_l2_transfer.hpp"
#include "moonolith_par_volume_surface_l2_transfer.hpp"
#include "par_moonolith_config.hpp"

#include <vector>

namespace utopia {

    namespace moonolith {

        // FIXME move somewhere public
        class FETransferOptions : public Configurable {
        public:
            virtual ~FETransferOptions() = default;

            void read(Input &in) override {
                if (!Options()  //
                         .add_option("verbose", verbose, "Verbose output during computation.")
                         .add_option("n_var",
                                     n_var,
                                     "Number of dimensions of vector function. Useful with tensor-product spaces.")
                         .add_option("use_reference_frame",
                                     use_reference_frame,
                                     "Use reference element frame for intersections.")
                         .add_option("clear_non_essential_matrices",
                                     clear_non_essential_matrices,
                                     "Keeps only the final transfer matrix in memory and deletes the rest.")
                         .add_option("export_tensors", export_tensors_, "Exports tensors to disk.")
                         .add_option("has_covering", has_covering, "Constructs lagrange multiplier in intersections.")
                         .parse(in)) {
                    return;
                }

                in.get("coupling", [this](Input &in) {
                    in.get_all([this](Input &in) {
                        int from = -1, to = -1;

                        in.get("from", from);
                        in.get("to", to);

                        assert(from != -1);
                        assert(to != -1);

                        tags.emplace_back(from, to);
                    });
                });
            }

            bool verbose{false};
            bool has_covering{true};
            int n_var{1};
            std::vector<std::pair<int, int>> tags;
            bool clear_non_essential_matrices{true};
            bool export_tensors_{false};
            bool use_reference_frame{false};
        };

        template <class Matrix>
        class FETransferData {
        public:
            /// From Petrov-Galerkin assembly
            std::shared_ptr<Matrix> coupling_matrix;

            /// Mass matrix assembly typically diagonal
            std::shared_ptr<Matrix> mass_matrix;

            /// Basis transformation used for second order fe
            std::shared_ptr<Matrix> basis_transform_matrix;

            /// Full transformation
            std::shared_ptr<Matrix> transfer_matrix;

            void create_matrices() {
                coupling_matrix = std::make_shared<Matrix>();
                mass_matrix = std::make_shared<Matrix>();
                basis_transform_matrix = std::make_shared<Matrix>();
                transfer_matrix = std::make_shared<Matrix>();
            }

            void clear_matrices() {
                coupling_matrix = nullptr;
                mass_matrix = nullptr;
                basis_transform_matrix = nullptr;
                transfer_matrix = nullptr;
            }

            inline bool empty() const { return !static_cast<bool>(transfer_matrix); }
        };

        class FETransferPrepareData {
        public:
            using Matrix_t = Traits<FunctionSpace>::Matrix;
            using Vector_t = Traits<FunctionSpace>::Vector;
            using Scalar_t = Traits<FunctionSpace>::Scalar;

            inline static void tensorize(const Matrix_t &T_x, const SizeType n_var, Matrix_t &T) {
                auto max_nnz = utopia::max_row_nnz(T_x);
                T.sparse(layout(T_x), max_nnz, max_nnz);

                assert(!empty(T));
                assert(T.row_range().extent() % n_var == 0);

                Write<Matrix_t> w(T);
                each_read(T_x, [&](const SizeType i, const SizeType j, const Scalar_t value) {
                    for (SizeType k = 0; k < n_var; ++k) {
                        T.set(i + k, j + k, value);
                    }
                });
            }

            template <class TransferAssembler>
            static bool apply(const FETransferOptions &opts, TransferAssembler &t, FETransferData<Matrix_t> &data) {
                UTOPIA_TRACE_REGION_BEGIN("FETransferPrepareData::apply");

                data.create_matrices();

                auto &B = *data.coupling_matrix;
                auto &D = *data.mass_matrix;
                auto &Q = *data.basis_transform_matrix;
                auto &T = *data.transfer_matrix;

                utopia::convert(t.buffers.B.get(), B);
                utopia::convert(t.buffers.D.get(), D);
                utopia::convert(t.buffers.Q.get(), Q);

                bool ok = false;
                if (!empty(Q)) {
                    m_utopia_warning_once("using sum(D, 1) instead of diag(D)");

                    // handle_constraints_pre_process(data);

                    Vector_t d_inv = sum(D, 1);

                    e_pseudo_inv(d_inv, d_inv, 1e-12);

                    Matrix_t D_tilde_inv = diag(d_inv);
                    Matrix_t T_x = Q * D_tilde_inv * B;

                    // handle_constraints_post_process(data, T_x);

                    if (opts.n_var == 1) {
                        T = T_x;
                        ok = true;
                        if (opts.clear_non_essential_matrices) {
                            data.coupling_matrix = nullptr;
                            data.mass_matrix = nullptr;
                            data.basis_transform_matrix = nullptr;
                        }

                    } else {
                        // assert(false && "IMPLEMENT ME");
                        tensorize(T_x, opts.n_var, T);
                    }
                }

                rename("transfer_matrix", *data.transfer_matrix);

                if (opts.export_tensors_) {
                    write("load_transfer_matrix.m", *data.transfer_matrix);
                }

                UTOPIA_TRACE_REGION_END("FETransferPrepareData::apply");
                return ok;
            }
        };

        template <int Dim, int DimFrom, int DimTo = Dim>
        class FETransferAux final {
        public:
            using Matrix_t = Traits<FunctionSpace>::Matrix;
            using Vector_t = Traits<FunctionSpace>::Vector;
            using Scalar_t = Traits<FunctionSpace>::Scalar;

            static bool apply(const FETransferOptions &opts,
                              const FunctionSpace &from,
                              const FunctionSpace &to,
                              FETransferData<Matrix_t> &data) {
                auto &m_from = *from.raw_type<Dim>();
                auto &m_to = *to.raw_type<Dim>();

                ::moonolith::Communicator comm = m_from.mesh().comm();
                comm.barrier();

                bool ok = true;

                if (opts.has_covering) {
                    ::moonolith::ParL2Transfer<Scalar_t,
                                               Dim,
                                               ::moonolith::StaticMax<DimFrom, 1>::value,
                                               ::moonolith::StaticMax<DimTo, 1>::value>
                        assembler(comm);

                    // assembler.use_reference_frame(opts.use_reference_frame);

                    if (opts.tags.empty()) {
                        if (!assembler.assemble(m_from, m_to)) {
                            ok = false;
                        }
                    } else {
                        if (!assembler.assemble(m_from, m_to, opts.tags)) {
                            ok = false;
                        }
                    }

                    if (ok) {
                        FETransferPrepareData::apply(opts, assembler, data);
                    }

                } else {
                    ::moonolith::ParL2TransferNoCovering<Scalar_t,
                                                         Dim,
                                                         ::moonolith::StaticMax<DimFrom, 1>::value,
                                                         ::moonolith::StaticMax<DimTo, 1>::value>
                        assembler(comm);

                    if (!assembler.assemble(m_from, m_to)) {
                        ok = false;
                    }

                    if (ok) {
                        FETransferPrepareData::apply(opts, assembler, data);
                    }
                }

                return ok;
            }
        };

        class InvalidFETransferAux {
        public:
            using Matrix_t = Traits<FunctionSpace>::Matrix;

            static bool apply(const FETransferOptions &opts,
                              const FunctionSpace &from,
                              const FunctionSpace &to,
                              FETransferData<Matrix_t> &data) {
                assert(false && "invalid dimensions");
                return false;
            }
        };

        template <int Dim>
        class LowerDimFETransferAux {
        public:
            using Matrix_t = Traits<FunctionSpace>::Matrix;
            using Vector_t = Traits<FunctionSpace>::Vector;
            using Scalar_t = Traits<FunctionSpace>::Scalar;

            static bool apply(const FETransferOptions &opts,
                              const FunctionSpace &from,
                              const FunctionSpace &to,
                              FETransferData<Matrix_t> &data) {
                ::moonolith::Communicator comm = from.mesh().comm().raw_comm();

                auto &m_from = *from.raw_type<Dim>();
                auto &m_to = *to.raw_type<Dim>();

                ::moonolith::ParVolumeSurfaceL2Transfer<Scalar_t, Dim> assembler(comm);

                if (opts.tags.empty()) {
                    if (!assembler.assemble(m_from, m_to)) {
                        return false;
                    }
                } else {
                    utopia::err() << "[Error] not implemented\n";
                    assert(false && "IMPLEMENT ME!!!");
                }

                FETransferPrepareData::apply(opts, assembler, data);
                return true;
            }
        };

        template <>
        class FETransferAux<1, 1, 0> final : public InvalidFETransferAux {};

        template <>
        class FETransferAux<0, 0, 0> final : public InvalidFETransferAux {};

        template <>
        class FETransferAux<2, 2, 1> final : public LowerDimFETransferAux<2> {};

        template <>
        class FETransferAux<3, 3, 2> final : public LowerDimFETransferAux<3> {};

        template <int Dim>
        class FETransferDispatch final {
        public:
            using Matrix_t = Traits<FunctionSpace>::Matrix;
            using Vector_t = Traits<FunctionSpace>::Vector;
            using Scalar_t = Traits<FunctionSpace>::Scalar;

            static bool apply(const FETransferOptions &opts,
                              const FunctionSpace &from,
                              const FunctionSpace &to,
                              FETransferData<Matrix_t> &data) {
                if (to.mesh().manifold_dimension() < Dim) {
                    if (from.mesh().manifold_dimension() < Dim) {
                        // surf to surf
                        assert(Dim > 1);
                        return FETransferAux<Dim, Dim - 1, Dim - 1>::apply(opts, from, to, data);
                    } else {
                        // Vol 2 surf
                        assert(Dim > 1);
                        return FETransferAux<Dim, Dim, Dim - 1>::apply(opts, from, to, data);
                    }
                } else {
                    // vol 2 vol
                    FETransferAux<Dim, Dim, Dim>::apply(opts, from, to, data);
                }

                return true;
            }
        };

        class FETransfer::Impl {
        public:
            FETransferOptions opts;
            FETransferData<Matrix> data;
        };

        FETransfer::FETransfer() : impl_(utopia::make_unique<Impl>()) {}
        FETransfer::~FETransfer() = default;

        void FETransfer::read(Input &in) { impl_->opts.read(in); }

        bool FETransfer::init(const std::shared_ptr<FunctionSpace> &from, const std::shared_ptr<FunctionSpace> &to) {
            UTOPIA_TRACE_REGION_BEGIN("FETransfer::init");

            clear();

            bool has_intersection = false;

            switch (from->mesh().spatial_dimension()) {
                case 1: {
                    has_intersection = FETransferDispatch<1>::apply(impl_->opts, *from, *to, impl_->data);
                    break;
                }
                case 2: {
                    has_intersection = FETransferDispatch<2>::apply(impl_->opts, *from, *to, impl_->data);
                    break;
                }
                case 3: {
                    has_intersection = FETransferDispatch<3>::apply(impl_->opts, *from, *to, impl_->data);
                    break;
                }
                default: {
                    assert(false && "Case not handled");
                    Utopia::Abort();
                }
            }

            UTOPIA_TRACE_REGION_END("FETransfer::init");
            return has_intersection;
        }

        void FETransfer::clear() { impl_->data.clear_matrices(); }

        bool FETransfer::empty() const { return impl_->data.empty(); }

        void FETransfer::describe(std::ostream &) const {}

        bool FETransfer::apply(const Vector &from, Vector &to) const {
            if (!empty()) {
                auto &op = *impl_->data.transfer_matrix;

                auto n_op = op.cols();
                auto n_from = from.size();

                auto n_var = n_from / n_op;

                assert(n_op * n_var == n_from);

                if (n_var == 1) {
                    to = (*impl_->data.transfer_matrix) * from;
                } else {
                    utopia::out() << "Tensorizing withing FETransfer::apply! n_var = " << n_var << "\n";

                    Vector scalar_from, scalar_to;
                    if (!utopia::empty(to)) {
                        to.set(0.0);
                    }

                    for (int c = 0; c < n_var; ++c) {
                        extract_component(from, n_var, c, scalar_from);
                        scalar_to = (*impl_->data.transfer_matrix) * scalar_from;
                        set_component(scalar_to, n_var, c, to);
                    }
                }

                return true;
            } else {
                return false;
            }
        }

        bool FETransfer::apply(const Matrix &to_matrix, Matrix &matrix_in_from_space) const {
            if (!empty()) {
                assert(!utopia::empty(*impl_->data.transfer_matrix));
                assert(!utopia::empty(to_matrix));
                matrix_in_from_space =
                    transpose(*impl_->data.transfer_matrix) * to_matrix * (*impl_->data.transfer_matrix);
                return true;
            } else {
                return false;
            }
        }

        bool FETransfer::apply_transpose(const Vector &from, Vector &to) const {
            if (!empty()) {
                to = transpose(*impl_->data.transfer_matrix) * from;
                return true;
            } else {
                return false;
            }
        }

        Size FETransfer::size() const {
            if (empty()) {
                return {0};
            } else {
                return impl_->data.coupling_matrix->size();
            }
        }

        Size FETransfer::local_size() const {
            if (empty()) {
                return {0};
            } else {
                return impl_->data.coupling_matrix->local_size();
            }
        }

        bool FETransfer::write(const Path &) const { return false; }

        FETransfer::Communicator &FETransfer::comm() {
            if (empty()) {
                assert(false);
                static Communicator self(Communicator::self());
                return self;
            } else {
                return impl_->data.coupling_matrix->comm();
            }
        }

        const FETransfer::Communicator &FETransfer::comm() const {
            if (empty()) {
                assert(false);
                static Communicator self(Communicator::self());
                return self;
            } else {
                return impl_->data.coupling_matrix->comm();
            }
        }

    }  // namespace moonolith

}  // namespace utopia

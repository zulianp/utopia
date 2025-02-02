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

            // constraint matrices (to be passed from outside)
            std::shared_ptr<Matrix> constraint_matrix_from;
            std::shared_ptr<Matrix> constraint_matrix_to;

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

            static void tensorize(const Matrix_t &T_x, const SizeType n_var, Matrix_t &T) {
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

            static void pre_constrain(FETransferData<Matrix_t> &data) {
                auto &B = *data.coupling_matrix;
                auto &D = *data.mass_matrix;

                if (data.constraint_matrix_to) {
                    auto &c_to = (*data.constraint_matrix_to);
                    D = transpose(c_to) * D * c_to;
                }

                if (data.constraint_matrix_from) {
                    auto &c_from = (*data.constraint_matrix_from);
                    if (data.constraint_matrix_to) {
                        // Both spaces are constrained
                        auto &c_to = (*data.constraint_matrix_to);
                        B = transpose(c_to) * B * c_from;
                    } else {
                        B = B * c_from;
                    }
                }
            }

            static void post_constrain(FETransferData<Matrix_t> &data) {
                // auto &Q = *data.basis_transform_matrix;
                auto &T = *data.transfer_matrix;

                if (data.constraint_matrix_to) {
                    auto &c_to = (*data.constraint_matrix_to);
                    Vector_t sum_c_to = sum(c_to, 1);

                    {
                        auto view = local_view_device(sum_c_to);
                        parallel_for(
                            local_range_device(sum_c_to), UTOPIA_LAMBDA(Size_t i) {
                                if (view.get(i) > 0.99) {
                                    view.set(i, 0.);
                                } else {
                                    view.set(i, 1.);
                                }
                            });
                    }

                    Matrix_t temp = diag(sum_c_to) * T;
                    T += temp + c_to * T;
                }
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

                    pre_constrain(data);

                    Vector_t d_inv = sum(D, 1);

                    e_pseudo_inv(d_inv, d_inv, 1e-12);

                    Matrix_t D_tilde_inv = diag(d_inv);
                    Matrix_t T_x = Q * D_tilde_inv * B;

                    // handle_constraints_post_process(data, T_x);

                    if (opts.chop_tol != 0.0) {
                        chop(T_x, opts.chop_tol);
                    }

                    post_constrain(data);

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
                              const FunctionSpace &from_and_to,
                              FETransferData<Matrix_t> &data) {
                auto &m_from_and_to = *from_and_to.raw_type<Dim>();
                ::moonolith::Communicator comm = m_from_and_to.mesh().comm();
                comm.barrier();

                ::moonolith::ParL2Transfer<Scalar_t, Dim, DimFrom, DimFrom> assembler(comm);
                bool ok = assembler.assemble(m_from_and_to, opts.tags, 1e-8);

                if (ok) {
                    FETransferPrepareData::apply(opts, assembler, data);
                }

                return ok;
            }

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

            static bool apply(const FETransferOptions &,
                              const FunctionSpace &,
                              const FunctionSpace &,
                              FETransferData<Matrix_t> &) {
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
                    // utopia::err() << "[Error] not implemented\n";
                    // assert(false && "IMPLEMENT ME!!!");

                    if (!assembler.assemble(m_from, m_to, opts.tags)) {
                        return false;
                    }
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

            static bool apply(const FETransferOptions &opts,
                              const FunctionSpace &from_and_to,
                              FETransferData<Matrix_t> &data) {
                if (from_and_to.mesh().manifold_dimension() < Dim) {
                    // surf to surf, avoid 0 for the moment (no point clouds)
                    static const int ManifoldDim = ::moonolith::StaticMax<Dim - 1, 1>::value;

                    assert(Dim > 1);
                    return FETransferAux<Dim, ManifoldDim, ManifoldDim>::apply(opts, from_and_to, data);
                } else {
                    // vol 2 vol
                    FETransferAux<Dim, Dim, Dim>::apply(opts, from_and_to, data);
                }

                return true;
            }
        };

        class FETransfer::Impl {
        public:
            FETransferOptions opts;
            FETransferData<Matrix> data;
        };

        void FETransfer::set_options(const FETransferOptions &options) {
            assert(impl_);
            if (impl_) {
                impl_->opts = options;
            }
        }

        void FETransfer::verbose(const bool val) {
            assert(impl_);
            if (impl_) {
                impl_->opts.verbose = val;
            }
        }

        FETransfer::FETransfer() : impl_(utopia::make_unique<Impl>()) {}
        FETransfer::~FETransfer() = default;

        void FETransfer::read(Input &in) { impl_->opts.read(in); }

        bool FETransfer::init(const std::shared_ptr<FunctionSpace> &from_and_to) {
            UTOPIA_TRACE_REGION_BEGIN("FETransfer::init");

            bool has_intersection = false;

            switch (from_and_to->mesh().spatial_dimension()) {
                case 1: {
                    has_intersection = FETransferDispatch<1>::apply(impl_->opts, *from_and_to, impl_->data);
                    break;
                }
                case 2: {
                    has_intersection = FETransferDispatch<2>::apply(impl_->opts, *from_and_to, impl_->data);
                    break;
                }
                case 3: {
                    has_intersection = FETransferDispatch<3>::apply(impl_->opts, *from_and_to, impl_->data);
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

        bool FETransfer::init(const std::shared_ptr<FunctionSpace> &from, const std::shared_ptr<FunctionSpace> &to) {
            UTOPIA_TRACE_REGION_BEGIN("FETransfer::init");

            ::moonolith::Moonolith::instance().verbose(Utopia::instance().verbose());

            assert(from->mesh().spatial_dimension() == to->mesh().spatial_dimension());

            if (from->mesh().spatial_dimension() != to->mesh().spatial_dimension()) {
                Utopia::Abort("Trying to transfer information between meshes with incompatibile embedding");
            }

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

            if (has_intersection) {
                // FIXME
                has_intersection = static_cast<bool>(impl_->data.transfer_matrix);
            }

            if (impl_->opts.verbose) {
                if (has_intersection) {
                    auto &T = *impl_->data.transfer_matrix;
                    auto &&comm = T.comm();

                    Scalar n_coupled_dofs = sum(T);
                    std::stringstream root_ss;
                    root_ss << "n_coupled_dofs in to_space:\t" << SizeType(n_coupled_dofs) << "\n";
                    root_ss << "global_nnz:\t" << T.global_nnz() << "\n";

                    if (impl_->opts.print_operator_imbalance) {
                        Matrix T_transpose = transpose(T);
                        std::vector<SizeType> balancing(comm.size(), 0);
                        balancing[comm.rank()] = T_transpose.local_nnz();

                        comm.sum(comm.size(), &balancing[0]);

                        if (comm.is_root()) {
                            root_ss << "Operator imbalance:\n";

                            for (int i = 0; i < comm.size(); ++i) {
                                if (balancing[i] > 0) {
                                    root_ss << "[" << i << " " << balancing[i] << "] ";
                                }
                            }

                            root_ss << "\n";
                        }
                    }

                    comm.root_print(root_ss.str(), utopia::out().stream());
                } else {
                    from->comm().root_print("NO INTERSECTION");
                }
            }

            UTOPIA_TRACE_REGION_END("FETransfer::init");

            if (impl_->opts.export_tensors) {
                ::utopia::write("load_transfer_matrix.m", *transfer_matrix());
            }

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

                auto n_var_from = n_from / n_op;

                assert(n_op * n_var_from == n_from);

                if (n_var_from == 1) {
                    to = (*impl_->data.transfer_matrix) * from;
                } else {
                    utopia::out() << "Tensorizing within FETransfer::apply! n_var_from = " << n_var_from << "\n";

                    Vector scalar_from, scalar_to;
                    if (!utopia::empty(to)) {
                        to.set(0.0);
                    }

                    SizeType n_to = to.size();

                    for (int c = 0; c < n_var_from; ++c) {
                        extract_component(from, n_var_from, c, scalar_from);
                        scalar_to = (*impl_->data.transfer_matrix) * scalar_from;
                        SizeType scalar_n_to = scalar_to.size();

                        if (scalar_n_to == n_to) {
                            copy_component(scalar_to, n_var_from, 0, c, to);
                        } else {
                            set_component(scalar_to, n_var_from, c, to);
                        }
                    }
                }

                return true;
            } else {
                return false;
            }
        }

        std::shared_ptr<FETransfer::Matrix> FETransfer::transfer_matrix() const { return impl_->data.transfer_matrix; }

        bool FETransfer::apply(const Matrix &to_matrix, Matrix &from_matrix, const bool reuse_matrix) const {
            if (!empty()) {
                assert(!utopia::empty(*impl_->data.transfer_matrix));
                assert(!utopia::empty(to_matrix));

                if (!utopia::empty(from_matrix) && reuse_matrix) {
                    ptap_reuse_matrix(to_matrix, *impl_->data.transfer_matrix, from_matrix);
                } else {
                    from_matrix = transpose(*impl_->data.transfer_matrix) * to_matrix * (*impl_->data.transfer_matrix);
                }
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
            if (empty() || !impl_->data.coupling_matrix) {
                static Communicator self(Communicator::self());
                return self;
            } else {
                return impl_->data.coupling_matrix->comm();
            }
        }

        const FETransfer::Communicator &FETransfer::comm() const {
            if (empty() || !impl_->data.coupling_matrix) {
                static Communicator self(Communicator::self());
                return self;
            } else {
                return impl_->data.coupling_matrix->comm();
            }
        }

        void FETransfer::set_constraint_matrix_from(const std::shared_ptr<Matrix> &constraint_matrix) {
            impl_->data.constraint_matrix_from = constraint_matrix;
        }

        void FETransfer::set_constraint_matrix_to(const std::shared_ptr<Matrix> &constraint_matrix) {
            impl_->data.constraint_matrix_to = constraint_matrix;
        }

        // void count_transpose_operator_imbalance(const Matrix &mat, IndexArray &array) {
        // auto &&comm = mat.comm();
        // array.resize(comm.size());

        //
        // }

        // void FETransfer::rebalance_the_from_system() {
        //     // TODO
        //     // Redistribute transfer matrix columns
        //     // Generate index set OR ranges for external routines
        //     // std::vector<long>
        // }

    }  // namespace moonolith

}  // namespace utopia

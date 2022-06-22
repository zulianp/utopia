#include "utopia_petsc_BDDOperator.hpp"

// Concrete includes
#include "utopia_DeviceView.hpp"
#include "utopia_TrivialPreconditioners.hpp"
#include "utopia_make_unique.hpp"

#include "utopia_petsc_CrsView.hpp"
#include "utopia_petsc_Factorization.hpp"
#include "utopia_petsc_Matrix.hpp"
#include "utopia_petsc_Matrix_impl.hpp"
#include "utopia_petsc_Vector.hpp"

#include "utopia_petsc_Utils.hpp"

#include <sstream>

namespace utopia {

    template <class Matrix, class Vector>
    class BDDOperator<Matrix, Vector>::Impl {
    public:
        static void parallel_to_serial(const Vector &x_from, Vector &x_to) {
            {
                // FIXME (copying stuff just because of abstractions!)
                auto x_from_view = local_view_device(x_from);
                auto x_to_view = local_view_device(x_to);
                parallel_for(
                    local_range_device(x_to),
                    UTOPIA_LAMBDA(const SizeType i) { x_to_view.set(i, x_from_view.get(i)); });
            }
        }

        static void serial_to_parallel(const Vector &x_from, Vector &x_to) { parallel_to_serial(x_from, x_to); }

        bool init(const std::shared_ptr<Matrix> &A_GG,
                  const std::shared_ptr<Matrix> &A_GI,
                  const std::shared_ptr<Matrix> &A_II,
                  const std::shared_ptr<Matrix> &A_IG) {
            A_GG_ = A_GG;
            A_GI_ = A_GI;
            A_II_ = A_II;
            A_IG_ = A_IG;

            UTOPIA_TRACE_REGION_BEGIN("BDDOperator::initialize::decomposition");
            // A_II_inv_ = std::make_shared<Factorization<Matrix, Vector>>("superlu", "lu");
            A_II_inv_ = std::make_shared<Factorization<Matrix, Vector>>("mumps", "cholesky");
            A_II_inv_->update(A_II_);
            UTOPIA_TRACE_REGION_END("BDDOperator::initialize::decomposition");
            return true;
        }

        bool init_rhs(const Vector &rhs_G, const Vector &rhs_I) {
            assert(A_II_inv_);

            // Compute rhs
            secant_G_ = std::make_shared<Vector>(layout(rhs_G));
            xL_ = std::make_shared<Vector>(row_layout(*A_GI_));
            rhsL_ = std::make_shared<Vector>(row_layout(*A_GI_));
            A_IG_x_ = std::make_shared<Vector>(row_layout(*A_II_));
            sol_I_ = std::make_shared<Vector>(row_layout(*A_II_));

            Vector inv_A_II_rhs_I(layout(rhs_I), 0);

            UTOPIA_TRACE_REGION_BEGIN("BDDOperator::init_rhs::apply");

            A_II_inv_->apply(rhs_I, inv_A_II_rhs_I);

            UTOPIA_TRACE_REGION_END("BDDOperator::init_rhs::apply");

            (*secant_G_) = rhs_G;

            Vector temp = (*A_GI_) * inv_A_II_rhs_I;

            assert(secant_G_->local_size() == temp.local_size());

            {
                auto sG_view = local_view_device(*secant_G_);
                auto temp_view = local_view_device(temp);

                parallel_for(
                    local_range_device(temp),
                    UTOPIA_LAMBDA(const SizeType i) { sG_view.set(i, sG_view.get(i) - temp_view.get(i)); });
            }

            return true;
        }

        bool apply(const Vector &x_G, Vector &rhs_G) const {
            parallel_to_serial(x_G, *xL_);

            *A_IG_x_ = (*A_IG_) * (*xL_);
            sol_I_->set(0);

            A_II_inv_->apply(*A_IG_x_, *sol_I_);

            (*rhsL_) = (*A_GI_) * (*sol_I_);

            if (empty(rhs_G)) {
                rhs_G.zeros(layout(x_G));
            }

            (*rhsL_) *= -1;

            serial_to_parallel(*rhsL_, rhs_G);
            rhs_G += (*A_GG_) * x_G;
            return true;
        }

        bool finalize(const Vector &x_G, const Vector &rhs_I, Vector &x_I) {
            parallel_to_serial(x_G, *xL_);

            *A_IG_x_ = (*A_IG_) * (*xL_);
            *rhsL_ = rhs_I;
            *rhsL_ -= *A_IG_x_;

            if (empty(x_I)) {
                x_I.zeros(layout(*rhsL_));
            } else {
                x_I.set(0);
            }

            A_II_inv_->apply(*rhsL_, x_I);
            return true;
        }

        void check(const Matrix &A) const {
            if (!debug) return;

            Matrix c_A = A;

            c_A.transform([](const Scalar &) { return 1; });

            Matrix diff = c_A - transpose(c_A);
            SizeType n_diff = norm1(diff);

            if (n_diff != 0) {
                std::stringstream ss;
                ss << "[Error] matrix does not have a symmetric graph. N diff = " << SizeType(n_diff) << "\n";
                A.comm().root_print(ss.str());
                A.comm().barrier();

                assert(false);

                Utopia::Abort("Quitting!");
            }
        }

        void select_constrained_rows(const Matrix &mat) {
            const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon();

            Vector mask(row_layout(mat), 0);

            if (selector.empty()) {
                selector.resize(mask.range().extent(), false);
            } else if (SizeType(selector.size()) != mask.range().extent()) {
                selector.resize(mask.range().extent());
                std::fill(std::begin(selector), std::end(selector), false);
            }

            {
                auto mask_view = view_device(mask);

                mat.read(UTOPIA_LAMBDA(const SizeType &i, const SizeType &j, const Scalar &value) {
                    if (i == j) return;

                    if (device::abs(value) > off_diag_tol) {
                        mask_view.atomic_add(i, 1);
                    }
                });
            }

            {
                auto mask_view = local_view_device(mask);

                parallel_for(
                    local_range_device(mask), UTOPIA_LAMBDA(const SizeType i) {
                        if (mask_view.get(i) == 0) {
                            selector[i] = true;
                        }
                    });

                if (verbose) {
                    parallel_for(
                        local_range_device(mask), UTOPIA_LAMBDA(const SizeType i) {
                            if (mask_view.get(i) == 0) {
                                mask_view.set(i, 1);
                            } else {
                                mask_view.set(i, 0);
                            }
                        });
                }
            }

            if (verbose) {
                SizeType n_linear_constraints = sum(mask);

                if (mat.comm().rank() == 0) {
                    utopia::out() << "n_linear_constraints: " << n_linear_constraints << "\n";
                }
            }
        }

        void fix_selected() {
            if (block_size <= 1) return;

            const SizeType n = selector.size();

            SizeType n_selected_before = 0;
            SizeType n_selected_after = 0;

            for (SizeType i = 0; i < n; i += block_size) {
                bool is_selected = false;

                for (int d = 0; d < block_size; ++d) {
                    bool is_selected_i = selector[i + d];
                    is_selected = is_selected || is_selected_i;

                    n_selected_before += is_selected_i;
                }

                if (is_selected) {
                    for (int d = 0; d < block_size; ++d) {
                        selector[i + d] = true;
                    }

                    n_selected_after += block_size;
                }
            }

            if (verbose) {
                std::stringstream ss;

                ss << "n_selected_before:\t" << n_selected_before << "\n";
                ss << "n_selected_after:\t" << n_selected_after << "\n";
                Communicator().root_print(ss.str());
            }
        }

        void init_interface(const Matrix &A) {
            UTOPIA_TRACE_REGION_BEGIN("BDDOperator::init_interface");

            check(A);
            if (handle_linear_constraints) {
                select_constrained_rows(A);
            }

            fix_selected();

            auto &&comm = A.comm();

            original_range = row_range(A);
            n_global = A.rows();

            SizeType n_local = original_range.extent();

            IndexSet min_idx(n_local, A.rows()), max_idx(n_local, 0);

            A.read([&](const SizeType i, const SizeType j, const Scalar) {
                min_idx[i - original_range.begin()] = std::min(SizeType(min_idx[i - original_range.begin()]), j);
                max_idx[i - original_range.begin()] = std::max(SizeType(max_idx[i - original_range.begin()]), j);
            });

            n_interface = 0;
            n_selected = 0;
            for (SizeType i = 0; i < n_local; ++i) {
                if (!original_range.inside(min_idx[i]) || !original_range.inside(max_idx[i])) {
                    ++n_interface;
                } else if (owned_is_selected(i)) {
                    ++n_selected;
                }
            }

            local_interface_idx.resize(n_interface + n_selected);
            std::fill(std::begin(local_interface_idx), std::end(local_interface_idx), 0);

            local_to_global.resize(n_interface + n_selected);
            std::fill(std::begin(local_to_global), std::end(local_to_global), 0);

            global_to_local.resize(n_local);
            std::fill(std::begin(global_to_local), std::end(global_to_local), -1);

            n_interface = 0;
            for (SizeType i = 0; i < n_local; ++i) {
                if (!original_range.inside(min_idx[i]) || !original_range.inside(max_idx[i])) {
                    local_to_global[n_interface] = (i + original_range.begin());
                    global_to_local[i] = n_interface;
                    local_interface_idx[n_interface++] = i;
                }
            }

            if (n_selected) {
                n_selected = 0;
                for (SizeType i = 0; i < n_local; ++i) {
                    if (global_to_local[i] == -1 && owned_is_selected(i)) {
                        local_to_global[n_interface + n_selected] = (i + original_range.begin());
                        global_to_local[i] = n_interface + n_selected;
                        local_interface_idx[n_interface + n_selected] = i;
                        ++n_selected;
                    }
                }
            }

            interface_offset = 0;
            selected_offset = 0;

            SizeType temp_buff_1[2] = {n_interface, n_selected};
            SizeType temp_buff_2[2] = {0, 0};
            comm.exscan_sum(temp_buff_1, temp_buff_2, 2);

            // comm.exscan_sum(&n_interface, &interface_offset, 1);

            interface_offset = temp_buff_2[0];
            selected_offset = temp_buff_2[1];

            UTOPIA_TRACE_REGION_END("BDDOperator::init_interface");
        }

        void build_A_II(const Matrix &A, Matrix &A_II) const {
            UTOPIA_TRACE_REGION_BEGIN("BDDOperator::build_A_II");

            Matrix A_II_view;
            local_block_view(A, A_II_view);

            A_II = A_II_view;
            set_zero_rows(A_II, local_interface_idx, 1);

            A_II.transform_ijv([&](const SizeType i, const SizeType j, const Scalar value) -> Scalar {
                bool is_removed = this->owned_is_in_G(j);
                if (is_removed) {
                    return i == j ? 1.0 : 0.;
                } else {
                    return value;
                }
            });

            UTOPIA_TRACE_REGION_END("BDDOperator::build_A_II");
        }

        void build_matrices(const Matrix &A, Matrix &A_GG, Matrix &A_IG, Matrix &A_GI) const {
            UTOPIA_TRACE_REGION_BEGIN("BDDOperator::build_matrices");

            auto &&comm = A.comm();

            auto r = row_range(A);
            SizeType n_local = r.extent();

            const PetscInt *ranges;
            MatGetOwnershipRanges(A.raw_type(), &ranges);

            SizeType n_col_dofs = A.cols();
            int comm_size = comm.size();
            auto find_rank = [=](const SizeType i) -> int {
                int rank = std::min(int(i * (float(comm_size) / n_col_dofs)), comm_size - 1);

                bool found = (i >= ranges[rank]) && (i < ranges[rank + 1]);

                while (!found) {
                    if (i < ranges[rank]) {
                        --rank;
                    } else if (i >= ranges[rank + 1]) {
                        ++rank;
                    } else {
                        assert(false);
                    }

                    assert(i >= 0);
                    assert(i < ranges[comm_size]);
                    assert(rank < comm_size);
                    assert(rank >= 0);

                    found = (i >= ranges[rank]) && (i < ranges[rank + 1]);
                }

                assert(found);

                return rank;
            };

            ////////////////////////////////////////////////////////////////////////

            std::vector<SizeType> counter(comm.size(), 0);
            std::vector<IndexSet> recv_list(comm.size()), send_list(comm.size());
            std::vector<SizeType> interface_ghosts;

            {
                Mat d, o;
                MatMPIAIJGetSeqAIJ(A.raw_type(), &d, &o, nullptr);

                utopia::PetscCrsView d_crs_view(d);
                utopia::PetscCrsView o_crs_view(o);

                PetscInt nghosts;
                const PetscInt *ghosts = nullptr;
                MatGetGhosts(A.raw_type(), &nghosts, &ghosts);
                interface_ghosts.resize(nghosts, -1);

                auto d_row_ptr = d_crs_view.row_ptr();
                auto d_colidx = d_crs_view.colidx();
                auto d_values = d_crs_view.values();

                auto o_row_ptr = o_crs_view.row_ptr();
                auto o_colidx = o_crs_view.colidx();
                auto o_values = o_crs_view.values();

                // Count incoming number of indices per process
                std::stringstream ss;
                for (SizeType i = 0; i < nghosts; ++i) {
                    int r = find_rank(ghosts[i]);
                    ++counter[r];
                }

                for (int r = 0; r < comm_size; ++r) {
                    if (r == comm.rank() || counter[r] == 0) continue;
                    recv_list[r].resize(counter[r]);
                }

                // Count outgoing number of indices per process
                std::fill(std::begin(counter), std::end(counter), 0);

                std::vector<SizeType> encountered_ranks;
                std::vector<bool> is_rank_counted(comm.size(), false);

                for (SizeType i = 0; i < n_interface; ++i) {
                    SizeType original_local_idx = local_interface_idx[i];
                    auto row_begin = o_row_ptr[original_local_idx];
                    auto row_end = o_row_ptr[original_local_idx + 1];

                    auto n_cols = row_end - row_begin;

                    if (SizeType(encountered_ranks.size()) < n_cols) {
                        encountered_ranks.resize(n_cols, -1);
                    }

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType col = o_colidx[k];
                        int rank = find_rank(ghosts[col]);

                        if (!is_rank_counted[rank]) {
                            ++counter[rank];
                        }

                        is_rank_counted[rank] = true;
                        encountered_ranks[k - row_begin] = rank;
                    }

                    for (SizeType k = 0; k < n_cols; ++k) {
                        assert(encountered_ranks[k] >= 0);
                        is_rank_counted[encountered_ranks[k]] = false;
                    }

#ifndef NDEBUG
                    for (auto irc : is_rank_counted) {
                        assert(!irc);
                    }
#endif
                }

                for (int r = 0; r < comm_size; ++r) {
                    if (r == comm.rank() || counter[r] == 0) continue;
                    send_list[r].reserve(counter[r]);
                }

                // Fill send_list with global indices
                for (SizeType i = 0; i < n_interface; ++i) {
                    SizeType original_local_idx = local_interface_idx[i];
                    SizeType new_global_idx = global_offset_G() + i;

                    auto row_begin = o_row_ptr[original_local_idx];
                    auto row_end = o_row_ptr[original_local_idx + 1];
                    auto n_cols = row_end - row_begin;

                    if (SizeType(encountered_ranks.size()) < n_cols) {
                        encountered_ranks.resize(n_cols, -1);
                    }

                    for (SizeType k = row_begin; k < row_end; ++k) {
                        SizeType col = o_colidx[k];
                        int rank = find_rank(ghosts[col]);

                        if (!is_rank_counted[rank]) {
                            send_list[rank].push_back(new_global_idx);
                        }

                        is_rank_counted[rank] = true;
                        encountered_ranks[k - row_begin] = rank;
                    }

                    for (SizeType k = 0; k < n_cols; ++k) {
                        assert(encountered_ranks[k] >= 0);
                        is_rank_counted[encountered_ranks[k]] = false;
                    }

#ifndef NDEBUG
                    for (auto irc : is_rank_counted) {
                        assert(!irc);
                    }
#endif
                }

                for (int r = 0; r < comm_size; ++r) {
                    if (send_list[r].empty()) {
                        assert(recv_list[r].empty());
                        continue;
                    }

                    int tag = 0;
                    MPI_Status status;
                    auto &send_buff = send_list[r];
                    auto &recv_buff = recv_list[r];

                    assert(!recv_buff.empty());

                    MPI_Sendrecv(send_buff.data(),
                                 send_buff.size(),
                                 MPIType<SizeType>::value(),
                                 r,
                                 tag,
                                 recv_buff.data(),
                                 recv_buff.size(),
                                 MPIType<SizeType>::value(),
                                 r,
                                 tag,
                                 comm.raw_comm(),
                                 &status);
                }

                SizeType linear_index = 0;
                for (int r = 0; r < comm_size; ++r) {
                    if (recv_list[r].empty()) continue;

                    for (auto idx : recv_list[r]) {
                        interface_ghosts[linear_index++] = idx;
                    }
                }

                ////////////////////////////////////////////////////////////////////////////////

                auto G_layout = layout(comm, local_size_G(), Traits::determine());
                auto GG_layout = square_matrix_layout(G_layout);

                IndexSet d_nnz(local_size_G(), 0), o_nnz(local_size_G(), 0);

                {
                    for (SizeType i = 0; i < local_size_G(); ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = d_row_ptr[original_local_idx];
                        auto row_end = d_row_ptr[original_local_idx + 1];

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];
                            if (this->owned_is_in_G(col)) {
                                ++d_nnz[i];
                            }
                        }
                    }

                    // count off proc nnz
                    for (SizeType i = 0; i < local_size_G(); ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = o_row_ptr[original_local_idx];
                        auto row_end = o_row_ptr[original_local_idx + 1];
                        auto n_cols = row_end - row_begin;
                        o_nnz[i] = n_cols;
                    }
                }

                ////////////////////////////////////////////////////////////////////////////////

                A_GG.sparse(GG_layout, d_nnz, o_nnz);

                {
                    Write<Matrix> w_A(A_GG);
                    for (SizeType i = 0; i < local_size_G(); ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = d_row_ptr[original_local_idx];
                        auto row_end = d_row_ptr[original_local_idx + 1];
                        SizeType new_row = global_offset_G() + i;

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];
                            SizeType new_col_local = global_to_local[col];

                            if (new_col_local != -1) {
                                SizeType new_col = global_offset_G() + new_col_local;
                                Scalar value = d_values[k];
                                A_GG.set(new_row, new_col, value);
                            }
                        }
                    }

                    // count off proc nnz
                    for (SizeType i = 0; i < local_size_G(); ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = o_row_ptr[original_local_idx];
                        auto row_end = o_row_ptr[original_local_idx + 1];
                        SizeType new_row = global_offset_G() + i;

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            Scalar value = o_values[k];
                            SizeType col = o_colidx[k];
                            SizeType new_col = interface_ghosts[col];
                            A_GG.set(new_row, new_col, value);
                        }
                    }
                }

                ////////////////////////////////////////////////////////////////////////////////

                auto GI_layout = serial_layout(local_size_G(), n_local);
                auto IG_layout = serial_layout(n_local, local_size_G());
                auto I_layout = serial_layout(n_local);

                std::fill(std::begin(d_nnz), std::end(d_nnz), 0);
                std::fill(std::begin(o_nnz), std::end(o_nnz), 0);

                {
                    for (SizeType i = 0; i < local_size_G(); ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = d_row_ptr[original_local_idx];
                        auto row_end = d_row_ptr[original_local_idx + 1];
                        // auto n_cols = row_end - row_begin;

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];
                            if (!this->owned_is_in_G(col)) {
                                ++d_nnz[i];
                            }
                        }
                    }

                    A_GI.sparse(GI_layout, d_nnz, o_nnz);

                    Write<Matrix> w(A_GI);

                    for (SizeType i = 0; i < local_size_G(); ++i) {
                        SizeType original_local_idx = local_interface_idx[i];
                        auto row_begin = d_row_ptr[original_local_idx];
                        auto row_end = d_row_ptr[original_local_idx + 1];

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];
                            if (!this->owned_is_in_G(col)) {
                                Scalar value = d_values[k];
                                A_GI.set(i, col, value);
                            }
                        }
                    }
                }

                /////////////////////////////////////////////////

                d_nnz.resize(n_local);
                o_nnz.resize(n_local);

                std::fill(std::begin(d_nnz), std::end(d_nnz), 0);
                std::fill(std::begin(o_nnz), std::end(o_nnz), 0);

                {
                    for (SizeType i = 0; i < n_local; ++i) {
                        auto row_begin = d_row_ptr[i];
                        auto row_end = d_row_ptr[i + 1];

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];

                            if (!this->owned_is_in_G(col)) {
                                ++d_nnz[i];
                            }
                        }
                    }

                    A_IG.sparse(IG_layout, d_nnz, o_nnz);
                    Write<Matrix> w_A(A_IG);

                    for (SizeType i = 0; i < n_local; ++i) {
                        if (this->owned_is_in_G(i)) continue;

                        auto row_begin = d_row_ptr[i];
                        auto row_end = d_row_ptr[i + 1];

                        for (SizeType k = row_begin; k < row_end; ++k) {
                            SizeType col = d_colidx[k];

                            SizeType new_col_local = global_to_local[col];

                            if (new_col_local != -1) {
                                Scalar value = d_values[k];
                                A_IG.set(i, new_col_local, value);
                            }
                        }
                    }
                }
            }

            UTOPIA_TRACE_REGION_END("BDDOperator::build_matrices");
        }

        inline bool owned_is_in_G(const SizeType i_owned) const {
            assert(i_owned < SizeType(global_to_local.size()));
            return (global_to_local[i_owned] != -1);
        }

        inline bool owned_is_selected(const SizeType i_owned) const {
            assert(selector.empty() || i_owned < SizeType(selector.size()));
            return (!selector.empty()) && selector[i_owned];
        }

        inline SizeType global_offset_G() const { return interface_offset + selected_offset; }
        inline SizeType local_size_G() const { return n_interface + n_selected; }

        std::shared_ptr<Matrix> A_GG_, A_GI_, A_II_, A_IG_;
        std::shared_ptr<Vector> secant_G_;
        std::shared_ptr<Factorization<Matrix, Vector>> A_II_inv_;
        std::shared_ptr<Vector> xL_, rhsL_, A_IG_x_, sol_I_;
        Vector b_I;

        // Interface support
        IndexSet local_interface_idx;
        IndexSet global_to_local;
        IndexSet local_to_global;
        SizeType interface_offset{0};
        SizeType selected_offset{0};
        SizeType n_interface{0};
        SizeType n_selected{0};
        Range original_range;
        SizeType n_global{0};

        Selector selector;

        std::string preconditioner_type{"inv_diag"};
        bool verbose{false};
        bool debug{false};
        bool handle_linear_constraints{true};

        int block_size{1};
    };

    template <class Matrix, class Vector>
    void BDDOperator<Matrix, Vector>::read(Input &in) {
        in.get("preconditioner_type", impl_->preconditioner_type);
        in.get("verbose", impl_->verbose);
        in.get("debug", impl_->debug);
        in.get("block_size", impl_->block_size);
        in.get("handle_linear_constraints", impl_->handle_linear_constraints);
    }

    template <class Matrix, class Vector>
    BDDOperator<Matrix, Vector>::BDDOperator() : impl_(utopia::make_unique<Impl>()) {}

    template <class Matrix, class Vector>
    BDDOperator<Matrix, Vector>::~BDDOperator() = default;

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::initialize(const std::shared_ptr<const Vector> &rhs) {
        UTOPIA_TRACE_REGION_BEGIN("BDDOperator::initialize(vector)");

        const auto &b = *rhs;
        auto &&comm = b.comm();

        const SizeType n_local = b.local_size();

        auto G_layout = layout(comm, impl_->local_size_G(), Traits::determine());
        auto I_layout = serial_layout(n_local);

        Vector b_G;
        Vector &b_I = impl_->b_I;

        {
            b_G.zeros(G_layout);

            auto b_G_view = local_view_device(b_G);
            auto b_view = local_view_device(b);

            for (SizeType i = 0; i < impl_->local_size_G(); ++i) {
                SizeType original_local_idx = impl_->local_interface_idx[i];
                b_G_view.set(i, b_view.get(original_local_idx));
            }

            b_I.zeros(I_layout);

            auto b_I_view = local_view_device(b_I);

            for (SizeType i = 0; i < n_local; ++i) {
                if (impl_->owned_is_in_G(i)) continue;
                b_I_view.set(i, b_view.get(i));
            }
        }

        bool ok = impl_->init_rhs(b_G, b_I);

        UTOPIA_TRACE_REGION_END("BDDOperator::initialize(vector)");
        return ok;
    }

    template <class Matrix, class Vector>
    typename BDDOperator<Matrix, Vector>::Selector &BDDOperator<Matrix, Vector>::selector() {
        return impl_->selector;
    }

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::initialize(const std::shared_ptr<const Matrix> &matrix) {
        UTOPIA_TRACE_REGION_BEGIN("BDDOperator::initialize(matrix)");

        auto A_GG_ptr = std::make_shared<Matrix>();
        auto A_GI_ptr = std::make_shared<Matrix>();
        auto A_IG_ptr = std::make_shared<Matrix>();
        auto A_II_ptr = std::make_shared<Matrix>();

        const Matrix &A = *matrix;
        Matrix &A_GG = *A_GG_ptr;
        Matrix &A_GI = *A_GI_ptr;
        Matrix &A_IG = *A_IG_ptr;
        Matrix &A_II = *A_II_ptr;

        impl_->init_interface(A);
        impl_->build_A_II(A, A_II);
        impl_->build_matrices(A, A_GG, A_IG, A_GI);

        bool ok = impl_->init(A_GG_ptr, A_GI_ptr, A_II_ptr, A_IG_ptr);

        UTOPIA_TRACE_REGION_END("BDDOperator::initialize(matrix)");
        return ok;
    }

    template <class Matrix, class Vector>
    const Vector &BDDOperator<Matrix, Vector>::righthand_side() const {
        return *impl_->secant_G_;
    }

    template <class Matrix, class Vector>
    std::shared_ptr<Matrix> BDDOperator<Matrix, Vector>::reduced_matrix() const {
        return impl_->A_GG_;
    }

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::apply(const Vector &x_G, Vector &rhs_G) const {
        return impl_->apply(x_G, rhs_G);
    }

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::finalize(const Vector &x_G, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BDDOperator::finalize");

        if (empty(x)) {
            x.zeros(layout(x_G.comm(), impl_->original_range.extent(), impl_->n_global));
        }

        bool ok = true;

        Vector x_I;
        if (!impl_->finalize(x_G, impl_->b_I, x_I)) {
            ok = false;
        } else {
            auto x_view = local_view_device(x);
            auto x_I_view = local_view_device(x_I);
            auto x_G_view = local_view_device(x_G);

            parallel_for(
                local_range_device(x), UTOPIA_LAMBDA(const SizeType i) { x_view.set(i, x_I_view.get(i)); });

            parallel_for(local_range_device(x_G), [&](const SizeType i) {
                SizeType i_local = impl_->local_to_global[i] - impl_->original_range.begin();
                x_view.set(i_local, x_G_view.get(i));
            });
        }

        UTOPIA_TRACE_REGION_END("BDDOperator::finalize");
        return ok;
    }

    template <class Matrix, class Vector>
    std::shared_ptr<Preconditioner<Vector>> BDDOperator<Matrix, Vector>::create_preconditioner() const {
        if (impl_->preconditioner_type == "inv") {
            if (impl_->verbose) comm().root_print("Using inv(A_GG) preconditioner");

            auto solver = std::make_shared<Factorization<Matrix, Vector>>();
            solver->update(impl_->A_GG_);
            return solver;
        } else if (impl_->preconditioner_type == "amg") {
            if (impl_->verbose) comm().root_print("Using amg(A_GG) preconditioner");

            auto solver = std::make_shared<KSPSolver<Matrix, Vector>>();
            solver->ksp_type("preonly");
            solver->pc_type("hypre");
            solver->update(impl_->A_GG_);
            solver->max_it(1);
            return solver;
        } else {
            if (impl_->verbose) comm().root_print("Using inv(diag(diag(A_GG)) preconditioner");

            auto d = std::make_shared<Vector>();
            this->diag(*d);
            *d = 1. / (*d);
            auto prec = std::make_shared<VectorPreconditioner<Vector>>(d);
            return prec;
        }
    }

    template <class Matrix, class Vector>
    typename BDDOperator<Matrix, Vector>::Layout BDDOperator<Matrix, Vector>::vector_layout() const {
        assert(impl_->A_GG_);
        return row_layout(*impl_->A_GG_);
    }

    template <class Matrix, class Vector>
    void BDDOperator<Matrix, Vector>::create_vector(Vector &x_G) {
        x_G.zeros(this->vector_layout());
    }

    template <class Matrix, class Vector>
    void BDDOperator<Matrix, Vector>::select(const Vector &x, Vector &x_G) const {
        auto &&comm = x.comm();
        auto G_layout = layout(comm, impl_->local_size_G(), Traits::determine());

        {
            x_G.zeros(G_layout);

            auto x_G_view = local_view_device(x_G);
            auto x_view = local_view_device(x);

            for (SizeType i = 0; i < impl_->local_size_G(); ++i) {
                SizeType original_local_idx = impl_->local_interface_idx[i];
                x_G_view.set(i, x_view.get(original_local_idx));
            }
        }
    }

    template <class Matrix, class Vector>
    void BDDOperator<Matrix, Vector>::diag(Vector &d) const {
        impl_->A_GG_->build_diag(d);
    }

    template <class Matrix, class Vector>
    Size BDDOperator<Matrix, Vector>::size() const {
        return impl_->A_GG_->size();
    }

    template <class Matrix, class Vector>
    Size BDDOperator<Matrix, Vector>::local_size() const {
        return impl_->A_GG_->local_size();
    }

    template <class Matrix, class Vector>
    typename Traits<Vector>::Communicator &BDDOperator<Matrix, Vector>::comm() {
        return impl_->A_GG_->comm();
    }

    template <class Matrix, class Vector>
    const typename Traits<Vector>::Communicator &BDDOperator<Matrix, Vector>::comm() const {
        return impl_->A_GG_->comm();
    }

    template class BDDOperator<PetscMatrix, PetscVector>;

}  // namespace utopia
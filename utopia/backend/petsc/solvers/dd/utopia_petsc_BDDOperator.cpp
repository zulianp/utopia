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

#include "utopia_petsc_Decompose.hpp"
#include "utopia_petsc_SchurComplement.hpp"

#include <sstream>

namespace utopia {

    template <class Matrix, class Vector>
    class BDDOperator<Matrix, Vector>::Impl {
    public:
        void init_skeleton_selection(const Matrix &mat) {
            SizeType local_rows = mat.local_rows();
            is_skeleton.resize(local_rows);

            if (!selector.empty()) {
                assert(SizeType(selector.size()) == local_rows);
                // Copy user selector as initial state
                is_skeleton = selector;
            } else {
                // Empty skeleton
                std::fill(std::begin(is_skeleton), std::end(is_skeleton), false);
            }

            // Add DD dofs to the skeleton
            auto rr = mat.row_range();

            IndexSet min_idx(local_rows, mat.rows());
            IndexSet max_idx(local_rows, 0);

            mat.read([&](const SizeType i, const SizeType j, const Scalar) {
                min_idx[i - rr.begin()] = std::min(SizeType(min_idx[i - rr.begin()]), j);
                max_idx[i - rr.begin()] = std::max(SizeType(max_idx[i - rr.begin()]), j);
            });

            for (SizeType i = 0; i < local_rows; ++i) {
                if (!rr.inside(min_idx[i]) || !rr.inside(max_idx[i])) {
                    is_skeleton[i] = true;
                }
            }
        }

        void extend_skeleton_selection_with_local_decomposition(const Matrix &mat) {
            if (num_blocks == 0) return;

#ifdef UTOPIA_ENABLE_METIS
            UTOPIA_TRACE_REGION_BEGIN("BDDOperator::extend_skeleton_with_local_blocks");

            Matrix local_block;
            local_block_view(mat, local_block);

            std::vector<int> partitions(local_block.rows(), -1);

            if (!decompose(local_block, num_blocks, &partitions[0])) {
                UTOPIA_TRACE_REGION_END("BDDOperator::extend_skeleton_with_local_blocks");
                Utopia::Abort("Failed to decompose into blocks");
            }

            auto mat_view = crs_view(local_block);

            int n_rows = mat_view.rows();
            int n_interface = 0;

            for (int r = 0; r < n_rows; ++r) {
                auto row = mat_view.row(r);

                for (int k = 0; k < int(row.length); ++k) {
                    auto c = row.colidx(k);

                    if (partitions[r] != partitions[c]) {
                        ++n_interface;
                        is_skeleton[r] = true;
                        break;
                    }
                }
            }

            if (verbose) {
                std::stringstream ss;
                ss << num_blocks << " blocks (" << n_interface << "/" << n_rows << ")\n";
                mat.comm().synched_print(ss.str(), utopia::out().stream());
            }

            UTOPIA_TRACE_REGION_END("BDDOperator::extend_skeleton_with_local_blocks");
#else
            UTOPIA_UNUSED(mat);
#endif  // UTOPIA_ENABLE_METIS
        }

        void initialize_dof_indices(const Matrix &mat) {
            const SizeType local_rows = mat.local_rows();

            SizeType n_skeleton = 0;
            for (SizeType i = 0; i < local_rows; ++i) {
                n_skeleton += is_skeleton[i];
            }

            const SizeType n_eliminated = local_rows - n_skeleton;

            skeleton_dofs.clear();
            eliminated_dofs.clear();

            skeleton_dofs.reserve(n_skeleton);
            eliminated_dofs.reserve(n_eliminated);

            auto rr = mat.row_range();
            for (SizeType i = 0; i < local_rows; ++i) {
                if (is_skeleton[i]) {
                    skeleton_dofs.push_back(i + rr.begin());
                }
            }

            for (SizeType i = 0; i < local_rows; ++i) {
                if (!is_skeleton[i]) {
                    eliminated_dofs.push_back(i + rr.begin());
                }
            }
        }

        void extend_skeleton_with_constrained_dofs(const Matrix &mat) {
            const Scalar off_diag_tol = std::numeric_limits<Scalar>::epsilon();

            Vector mask(row_layout(mat), 0);

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

                // FIXME ref capture
                parallel_for(local_range_device(mask), [&](const SizeType i) {
                    if (mask_view.get(i) == 0) {
                        is_skeleton[i] = true;
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
                std::stringstream ss;
                ss << "n_linear_constraints: " << n_linear_constraints << "\n";
                mat.comm().root_print(ss.str(), utopia::out().stream());
            }
        }

        void fix_skeleton_wrt_block_size(const Matrix &mat) {
            if (block_size <= 1) return;

            const SizeType n = is_skeleton.size();

            SizeType n_selected_before = 0;
            SizeType n_selected_after = 0;

            for (SizeType i = 0; i < n; i += block_size) {
                bool is_selected = false;

                for (int d = 0; d < block_size; ++d) {
                    bool is_selected_i = is_skeleton[i + d];
                    is_selected = is_selected || is_selected_i;

                    n_selected_before += is_selected_i;
                }

                if (is_selected) {
                    for (int d = 0; d < block_size; ++d) {
                        is_skeleton[i + d] = true;
                    }

                    n_selected_after += block_size;
                }
            }

            if (verbose) {
                auto &&comm = mat.comm();

                n_selected_before = comm.sum(n_selected_before);
                n_selected_after = comm.sum(n_selected_after);

                std::stringstream ss;
                ss << "Skeleton (before -> after): " << n_selected_before << " -> " << n_selected_after << " / "
                   << mat.rows() << " dofs\n";

                comm.root_print(ss.str());
            }
        }

        bool initialize_operator(const Matrix &mat) {
            init_skeleton_selection(mat);
            extend_skeleton_with_constrained_dofs(mat);
            extend_skeleton_selection_with_local_decomposition(mat);
            fix_skeleton_wrt_block_size(mat);
            initialize_dof_indices(mat);

            return schur_complement.initialize_from_selection(mat, eliminated_dofs);
        }

        bool initialize_righthand_side(const Vector &rhs) {
            return schur_complement.apply_righthand_side(rhs, rhs_eliminated, rhs_skeleton);
        }

        bool apply(const Vector &x_skeleton, Vector &y_skeleton) const {
            return schur_complement.apply(x_skeleton, y_skeleton);
        }

        bool finalize(const Vector &x_skeleton, Vector &x) {
            return schur_complement.finalize_from_eliminated_rhs(rhs_eliminated, x_skeleton, x);
        }

        Vector rhs_eliminated;
        Vector rhs_skeleton;

        IndexSet skeleton_dofs;
        IndexSet eliminated_dofs;

        // User selector
        Selector selector;

        // Internal selector
        Selector is_skeleton;

        bool verbose{false};
        // bool debug{false};
        // bool handle_linear_constraints{true};

        int block_size{1};
        int num_blocks{0};

        std::string preconditioner_type{"inv_diag"};

        SchurComplement<Matrix> schur_complement;
    };

    template <class Matrix, class Vector>
    void BDDOperator<Matrix, Vector>::read(Input &in) {
        in.get("preconditioner_type", impl_->preconditioner_type);
        in.get("verbose", impl_->verbose);
        // in.get("debug", impl_->debug);
        in.get("block_size", impl_->block_size);
        // in.get("handle_linear_constraints", impl_->handle_linear_constraints);
        in.get("num_blocks", impl_->num_blocks);
    }

    template <class Matrix, class Vector>
    BDDOperator<Matrix, Vector>::BDDOperator() : impl_(utopia::make_unique<Impl>()) {}

    template <class Matrix, class Vector>
    BDDOperator<Matrix, Vector>::~BDDOperator() = default;

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::initialize(const std::shared_ptr<const Vector> &rhs) {
        UTOPIA_TRACE_REGION_BEGIN("BDDOperator::initialize(vector)");
        assert(rhs);
        bool ok = impl_->initialize_righthand_side(*rhs);
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
        bool ok = impl_->initialize_operator(*matrix);
        UTOPIA_TRACE_REGION_END("BDDOperator::initialize(matrix)");
        return ok;
    }

    template <class Matrix, class Vector>
    const Vector &BDDOperator<Matrix, Vector>::righthand_side() const {
        return impl_->rhs_skeleton;
    }

    template <class Matrix, class Vector>
    std::shared_ptr<Matrix> BDDOperator<Matrix, Vector>::reduced_matrix() const {
        return impl_->schur_complement.reduced_matrix();
    }

    template <class Matrix, class Vector>
    const typename BDDOperator<Matrix, Vector>::IndexSet &BDDOperator<Matrix, Vector>::skeleton_dofs() const {
        return impl_->skeleton_dofs;
    }

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::apply(const Vector &x_skeleton, Vector &rhs_skeleton) const {
        UTOPIA_TRACE_REGION_BEGIN("BDDOperator::apply");

        bool ok = impl_->apply(x_skeleton, rhs_skeleton);

        UTOPIA_TRACE_REGION_END("BDDOperator::apply");
        return ok;
    }

    template <class Matrix, class Vector>
    bool BDDOperator<Matrix, Vector>::finalize(const Vector &x_skeleton, Vector &x) {
        UTOPIA_TRACE_REGION_BEGIN("BDDOperator::finalize");

        bool ok = impl_->finalize(x_skeleton, x);

        UTOPIA_TRACE_REGION_END("BDDOperator::finalize");
        return ok;
    }

    template <class Matrix, class Vector>
    std::shared_ptr<Preconditioner<Vector>> BDDOperator<Matrix, Vector>::create_preconditioner() const {
        if (impl_->preconditioner_type == "inv") {
            if (impl_->verbose) comm().root_print("Using inv(A_GG) preconditioner");

            auto solver = std::make_shared<Factorization<Matrix, Vector>>();
            solver->update(reduced_matrix());
            return solver;
        } else if (impl_->preconditioner_type == "amg") {
            if (impl_->verbose) comm().root_print("Using amg(A_GG) preconditioner");

            auto solver = std::make_shared<KSPSolver<Matrix, Vector>>();
            solver->ksp_type(KSPPREONLY);
            solver->pc_type(PCHYPRE);
            solver->update(reduced_matrix());
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
        return row_layout(*impl_->schur_complement.reduced_matrix());
    }

    template <class Matrix, class Vector>
    void BDDOperator<Matrix, Vector>::create_vector(Vector &x_skeleton) {
        x_skeleton.zeros(this->vector_layout());
    }

    template <class Matrix, class Vector>
    void BDDOperator<Matrix, Vector>::select(const Vector &x, Vector &x_skeleton) const {
        impl_->schur_complement.select(x, x_skeleton);
    }

    template <class Matrix, class Vector>
    void BDDOperator<Matrix, Vector>::diag(Vector &d) const {
        impl_->schur_complement.reduced_matrix()->build_diag(d);
    }

    template <class Matrix, class Vector>
    Size BDDOperator<Matrix, Vector>::size() const {
        return impl_->schur_complement.size();
    }

    template <class Matrix, class Vector>
    Size BDDOperator<Matrix, Vector>::local_size() const {
        return impl_->schur_complement.local_size();
    }

    template <class Matrix, class Vector>
    typename Traits<Vector>::Communicator &BDDOperator<Matrix, Vector>::comm() {
        return impl_->schur_complement.comm();
    }

    template <class Matrix, class Vector>
    const typename Traits<Vector>::Communicator &BDDOperator<Matrix, Vector>::comm() const {
        return impl_->schur_complement.comm();
    }

    template class BDDOperator<PetscMatrix, PetscVector>;

}  // namespace utopia

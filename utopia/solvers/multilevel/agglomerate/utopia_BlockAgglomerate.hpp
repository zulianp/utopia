#ifndef UTOPIA_BLOCK_AGGLOMERATE_HPP
#define UTOPIA_BLOCK_AGGLOMERATE_HPP

#include "utopia_Base.hpp"

#ifdef UTOPIA_WITH_PETSC

#include "utopia_AlgebraicMultigrid.hpp"

#include "utopia_CRSMatrix.hpp"
#include "utopia_petsc_ILUDecompose.hpp"

namespace utopia {

    template <class Matrix, int BlockSize, int Backend = Traits<Matrix>::Backend>
    class BlockAgglomerate final : public MatrixAgglomerator<Matrix> {
    public:
        using Traits = utopia::Traits<Matrix>;
        using Comm = typename Traits::Communicator;
        using Scalar = typename Traits::Scalar;
        using SizeType = typename Traits::SizeType;
        using IndexArray = typename Traits::IndexArray;
        using ScalarArray = typename Traits::ScalarArray;
        using Vector = typename Traits::Vector;
        using Transfer = utopia::Transfer<Matrix, Vector>;

        using BlockMatrix = utopia::CRSMatrix<std::vector<Scalar>, std::vector<SizeType>, BlockSize>;

        BlockAgglomerate *clone() const override { return new BlockAgglomerate(*this); }

        void read(Input &in) override {
            in.get("bmax", bmax_);
            in.get("weight", weight_);
            in.get("verbose", verbose_);
            in.get("component", component_);
        }

        inline Scalar block_weight(const Scalar *block) const { return block[component_ * BlockSize + component_]; }

        std::shared_ptr<Transfer> create_transfer(const Matrix &in) override {
            UTOPIA_TRACE_REGION_BEGIN("BlockAgglomerate::create_prolongator");
            using namespace utopia;

            auto prolongator = std::make_shared<Matrix>();

            auto rr = row_range(in);

            Matrix in_local;
            local_block_view(in, in_local);

            BlockMatrix block_mat;
            crs_block_matrix(in_local, block_mat);

            SizeType n_blocks = block_mat.rows();

            IndexArray parent(n_blocks, -1);
            SizeType n_coarse_rows = 0;

            ScalarArray a_max(n_blocks, 0);

            for (SizeType block_i = 0; block_i < n_blocks; ++block_i) {
                auto row_view = block_mat.row(block_i);

                for (SizeType k = 0; k < row_view.n_blocks(); ++k) {
                    if (row_view.colidx(k) == block_i) continue;
                    a_max[block_i] = device::max(a_max[block_i], device::abs(block_weight(row_view.block(k))));
                }
            }

            {
                for (SizeType block_i = 0; block_i < n_blocks; ++block_i) {
                    if (parent[block_i] != -1) continue;
                    parent[block_i] = n_coarse_rows;

                    auto row_view = block_mat.row(block_i);
                    SizeType n_values = row_view.n_blocks();

                    Scalar a_max_i = a_max[block_i];

                    SizeType count = 0;
                    for (SizeType k = 0; k < n_values; ++k) {
                        const SizeType block_j = row_view.colidx(k);
                        const Scalar a_ij = block_weight(row_view.block(k));

                        if (weight_ * a_max_i < device::abs(a_ij) && parent[block_j] == -1) {
                            parent[block_j] = n_coarse_rows;
                            count++;
                        }

                        if (count >= bmax_) {
                            break;
                        }
                    }

                    if (count == 0) {
                        // This node is not a cluster center
                        parent[block_i] = -1;

                        Scalar max_aij = 0;
                        SizeType arg_max_j = -1;
                        for (SizeType k = 0; k < n_values; ++k) {
                            const SizeType block_j = row_view.colidx(k);
                            const Scalar a_ij = std::abs(block_weight(row_view.block(k)));

                            if (a_ij > max_aij) {
                                arg_max_j = block_j;
                                max_aij = a_ij;
                            }
                        }

                        if (arg_max_j != -1) {
                            parent[block_i] = parent[arg_max_j];
                        }

                    } else {
                        ++n_coarse_rows;
                    }
                }

                SizeType n_not_aggr = 0;
                for (SizeType block_i = 0; block_i < n_blocks; ++block_i) {
                    if (parent[block_i] != -1) continue;

                    parent[block_i] = n_coarse_rows++;
                    ++n_not_aggr;
                }

                if (verbose_) {
                    in.comm().synched_print("n_not_aggr: " + std::to_string(n_not_aggr) +
                                            ", n_coarse_rows: " + std::to_string(n_coarse_rows) + "/" +
                                            std::to_string(in.local_rows()) + "\n");
                }
            }

            auto pl = layout(in.comm(), in.local_rows(), n_coarse_rows * BlockSize, in.rows(), Traits::determine());
            prolongator->sparse(pl, 1, 1);

            {
                Write<Matrix> w(*prolongator);
                auto coarse_offset = prolongator->col_range().begin();

                for (SizeType block_i = 0; block_i < n_blocks; ++block_i) {
                    for (int sub_i = 0; sub_i < BlockSize; ++sub_i) {
                        prolongator->set(rr.begin() + block_i * BlockSize + sub_i,
                                         coarse_offset + parent[block_i] * BlockSize + sub_i,
                                         1.0);
                    }
                }
            }

            UTOPIA_TRACE_REGION_END("BlockAgglomerate::create_prolongator");
            return std::make_shared<IPRTransfer<Matrix, Vector>>(prolongator);
        }

        inline void verbose(const bool val) { verbose_ = val; }

    private:
        SizeType bmax_{3};
        // SizeType bmin_{5};
        SizeType component_{BlockSize - 1};
        Scalar weight_{1. / 3};
        bool verbose_{false};
    };

}  // namespace utopia

#endif
#endif  // UTOPIA_BLOCK_AGGLOMERATE_HPP

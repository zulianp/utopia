#ifndef UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_HPP
#define UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_HPP

#include "utopia_ILU.hpp"
#include "utopia_MonotoneMultigrid.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_BoundingBoxMultiLevel.hpp"

#ifdef UTOPIA_WITH_BLAS
#include "utopia_RASPatchSmoother.hpp"
#include "utopia_blas.hpp"
#include "utopia_blas_Array.hpp"
#endif  // UTOPIA_WITH_BLAS

#include <memory>
#include <vector>

namespace utopia {

    template <class FunctionSpace>
    class MonotoneSemiGeometricMultigrid
        : public QPSolver<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector> {
    public:
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using IndexSet = typename Traits<FunctionSpace>::IndexSet;
        using Comm = typename Traits<FunctionSpace>::Communicator;

        using Super = utopia::QPSolver<Matrix, Vector>;
        using IPTruncatedTransfer = utopia::IPTruncatedTransfer<Matrix, Vector>;
        using IPRTransfer = utopia::IPRTransfer<Matrix, Vector>;
        using IPTransfer = utopia::IPTransfer<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using MonotoneMultigrid = utopia::MonotoneMultigrid<Matrix, Vector>;
        using Smoother = utopia::IterativeSolver<Matrix, Vector>;
        using QPSmoother = utopia::QPSolver<Matrix, Vector>;

        using CoarseSpaceGen = utopia::BoundingBoxMultiLevelFunctionSpaceGenerator<FunctionSpace>;

        void generate_coarse_spaces(Input &in, FunctionSpace &finest) {
            int n_auto_levels = space_hierarchy_.empty() ? 1 : 0;
            in.get("n_auto_levels", n_auto_levels);
            if (n_auto_levels == 0) return;

            std::vector<std::shared_ptr<FunctionSpace>> coarse_spaces;
            CoarseSpaceGen gen;
            gen.read(in);

            if (gen.generate_spaces(finest, n_auto_levels, coarse_spaces)) {
                space_hierarchy_.insert(space_hierarchy_.begin(), coarse_spaces.begin(), coarse_spaces.end());
            } else {
                Utopia::Abort("Failed to create coarse spaces!");
            }
        }

        void read(Input &in) override {
            Super::read(in);

            if (fine_space_) {
                block_size_ = fine_space_->n_var();
            }

            bool export_coarse_meshes = false;
            in.get("export_coarse_meshes", export_coarse_meshes);
            in.get("block_size", block_size_);
            in.get("use_patch_smoother", use_patch_smoother_);

            in.get("coarse_spaces", [this](Input &array_node) {
                array_node.get_all([this](Input &node) {
                    auto space = std::make_shared<FunctionSpace>();
                    space->read(node);

                    assert(!space->empty());
                    space_hierarchy_.push_back(space);
                });
            });

            in.get("clear_spaces_after_init", clear_spaces_after_init_);

            if (space_hierarchy_.empty()) {
                generate_coarse_spaces(in, *fine_space_);
            } else {
                generate_coarse_spaces(in, *space_hierarchy_[0]);
            }

            if (export_coarse_meshes) {
                int num_space = 0;
                for (auto &s : space_hierarchy_) {
                    s->mesh().write("mesh_" + std::to_string(num_space++) + ".e");
                }
            }

            const int n_levels = space_hierarchy_.size();

            if (n_levels > 0) {
                make_algo();
                algorithm_->read(in);
            }
        }

        MonotoneSemiGeometricMultigrid *clone() const override {
            auto ptr = utopia::make_unique<MonotoneSemiGeometricMultigrid>();

            if (this->algorithm_) {
                ptr->algorithm_ = std::unique_ptr<MonotoneMultigrid>(this->algorithm_->clone());
            }

            return ptr.release();
        }

        MonotoneSemiGeometricMultigrid() : algorithm_() {}

        ~MonotoneSemiGeometricMultigrid() override = default;

        /*! @brief if overriden the subclass has to also call this one first
         */
        void update(const std::shared_ptr<const Matrix> &op) override {
            assert(algorithm_);
            assert(fine_space_);

            if (!is_initialized_) {
                assert(!space_hierarchy_.empty());
                init();
                is_initialized_ = true;
            }

            algorithm_->update(op);
        }

        bool init() {
            const int n_levels = space_hierarchy_.size();

            if (n_levels == 0) {
                assert(false);
                return false;
            }

            if (!fine_space_) {
                assert(false);
                Utopia::Abort("MonotoneSemiGeometricMultigrid fine space must be set!");
            }

            if (empty()) {
                make_algo();
            }

            InputParameters params;
            params.set("n_var", fine_space_->n_var());
            params.set("chop_tol", 1e-8);

            std::vector<std::shared_ptr<Transfer>> transfers;
            FETransfer<FunctionSpace> transfer;
            transfer.read(params);

            for (int l = 0; l < n_levels - 1; ++l) {
                if (!transfer.init(space_hierarchy_[l], space_hierarchy_[l + 1])) {
                    assert(false);
                    Utopia::Abort("MonotoneSemiGeometricMultigrid failed to set-up transfer operator!");
                }

                // transfers.push_back(transfer.template build_transfer<IPRTransfer>());
                transfers.push_back(transfer.template build_transfer<IPTransfer>());
            }

            if (!transfer.init(space_hierarchy_[n_levels - 1], fine_space_)) {
                assert(false);
                Utopia::Abort("MonotoneSemiGeometricMultigrid failed to set-up transfer operator!");
            }

            auto t = transfer.template build_transfer<IPTruncatedTransfer>();
            fine_space_->apply_constraints(*t->I_ptr(), 0.0);

            // rename("T", const_cast<Matrix &>(t->I()));
            // write("load_T.m", t->I());

            transfers.push_back(t);
            algorithm_->set_transfer_operators(transfers);

            if (clear_spaces_after_init_) {
                space_hierarchy_.clear();
            }

            return true;
        }

        bool apply(const Vector &rhs, Vector &x) override {
            assert(algorithm_);
            if (!algorithm_) return false;

            algorithm_->set_box_constraints(this->get_box_constraints());
            return algorithm_->apply(rhs, x);
        }

        inline void set_fine_space(const std::shared_ptr<FunctionSpace> &fine_space) { fine_space_ = fine_space; }
        inline bool empty() const { return !algorithm_; }

    private:
        std::unique_ptr<MonotoneMultigrid> algorithm_;

        std::shared_ptr<FunctionSpace> fine_space_;
        std::vector<std::shared_ptr<FunctionSpace>> space_hierarchy_;
        bool is_initialized_{false};
        bool clear_spaces_after_init_{false};
        bool use_patch_smoother_{false};
        int block_size_{1};

        void make_algo() {
            InputParameters params;
            params.set("block_size", block_size_);

            std::shared_ptr<QPSmoother> fine_smoother;
            std::shared_ptr<Smoother> coarse_smoother;

#ifdef UTOPIA_WITH_BLAS
            if (use_patch_smoother_) {
                fine_smoother = std::make_shared<RASPatchSmoother<Matrix, utopia::BlasMatrixd>>();
                coarse_smoother = std::make_shared<RASPatchSmoother<Matrix, utopia::BlasMatrixd>>();
            } else
#endif  // UTOPIA_WITH_BLAS
            {
                fine_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
                coarse_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            }

            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>();

            fine_smoother->read(params);
            coarse_smoother->read(params);
            algorithm_ = utopia::make_unique<MonotoneMultigrid>(fine_smoother, coarse_smoother, direct_solver);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_HPP

#ifndef UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_NEW_HPP
#define UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_NEW_HPP

#include "utopia_ILU.hpp"
#include "utopia_Multigrid.hpp"
#include "utopia_ProjectedGaussSeidelNew.hpp"

#include "utopia_fe_Core.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class FunctionSpace>
    class SemiGeometricMultigridNew
        : public IterativeSolver<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector> {
    public:
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using IndexSet = typename Traits<FunctionSpace>::IndexSet;
        using Comm = typename Traits<FunctionSpace>::Communicator;

        using Super = utopia::IterativeSolver<Matrix, Vector>;
        using IPRTruncatedTransfer = utopia::IPRTruncatedTransfer<Matrix, Vector>;
        using IPRTransfer = utopia::IPRTransfer<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using Multigrid = utopia::Multigrid<Matrix, Vector>;

        void read(Input &in) override {
            Super::read(in);

            bool export_coarse_meshes = false;
            in.get("export_coarse_meshes", export_coarse_meshes);

            in.get("coarse_spaces", [this, export_coarse_meshes](Input &array_node) {
                array_node.get_all([this, export_coarse_meshes](Input &node) {
                    auto space = std::make_shared<FunctionSpace>();
                    space->read(node);

                    if (export_coarse_meshes) {
                        space->mesh().write("mesh_" + std::to_string(space_hierarchy_.size()) + ".e");
                    }

                    assert(!space->empty());
                    space_hierarchy_.push_back(space);
                });
            });

            const int n_levels = space_hierarchy_.size();

            if (n_levels > 0) {
                make_algo();
                algorithm_->read(in);
            }
        }

        SemiGeometricMultigridNew *clone() const override {
            auto ptr = utopia::make_unique<SemiGeometricMultigridNew>();

            if (this->algorithm_) {
                ptr->algorithm_ = std::unique_ptr<Multigrid>(this->algorithm_->clone());
            }

            return ptr.release();
        }

        SemiGeometricMultigridNew() : algorithm_() {}

        ~SemiGeometricMultigridNew() override = default;

        /*! @brief if overriden the subclass has to also call this one first
         */
        void update(const std::shared_ptr<const Matrix> &op) override {
            assert(algorithm_);
            assert(!space_hierarchy_.empty());
            assert(fine_space_);

            if (!is_initialized_) {
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
                Utopia::Abort("SemiGeometricMultigridNew fine space must be set!");
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
                    Utopia::Abort("SemiGeometricMultigridNew failed to set-up transfer operator!");
                }

                transfers.push_back(transfer.template build_transfer<IPRTransfer>());
            }

            if (!transfer.init(space_hierarchy_[n_levels - 1], fine_space_)) {
                assert(false);
                Utopia::Abort("SemiGeometricMultigridNew failed to set-up transfer operator!");
            }

            auto t = transfer.template build_transfer<IPRTruncatedTransfer>();

            // rename("T", const_cast<Matrix &>(t->I()));
            // write("load_T.m", t->I());

            transfers.push_back(t);
            algorithm_->set_transfer_operators(transfers);
            return true;
        }

        bool apply(const Vector &rhs, Vector &x) override {
            assert(algorithm_);
            if (!algorithm_) return false;
            return algorithm_->apply(rhs, x);
        }

        inline void set_fine_space(const std::shared_ptr<FunctionSpace> &fine_space) { fine_space_ = fine_space; }
        inline bool empty() const { return !algorithm_; }

    private:
        std::unique_ptr<Multigrid> algorithm_;

        std::shared_ptr<FunctionSpace> fine_space_;
        std::vector<std::shared_ptr<FunctionSpace>> space_hierarchy_;
        bool is_initialized_{false};

        void make_algo() {
            InputParameters params;
            params.set("block_size", fine_space_->n_var());

            auto smoother = std::make_shared<SOR<Matrix, Vector>>();
            auto direct_solver = std::make_shared<KSPSolver<Matrix, Vector>>();

            smoother->read(params);
            algorithm_ = utopia::make_unique<Multigrid>(smoother, direct_solver);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_NEW_HPP
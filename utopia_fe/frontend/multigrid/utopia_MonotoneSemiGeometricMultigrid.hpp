#ifndef UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_HPP
#define UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_HPP

#include "utopia_ILU.hpp"
#include "utopia_MonotoneMultigrid.hpp"
#include "utopia_SolverType.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_BoundingBoxMultiLevel.hpp"
#include "utopia_SemiGeometricHierarchy.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class FunctionSpace>
    class MonotoneSemiGeometricMultigrid
        : public QPSolver<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector>,
          public SemiGeometricHierarchy<FunctionSpace> {
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

        inline static constexpr SolverType solver_type() { return "MonotoneSemiGeometricMultigrid"; }
        inline static constexpr SolverType backend() { return Solver::any_backend(); }

        void read(Input &in) override {
            Super::read(in);

            this->read_hierarchy(in);

            const int n_levels = this->n_coarse_spaces();

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

            if (!is_initialized_) {
                init();
                is_initialized_ = true;
            }

            algorithm_->update(op);
        }

        bool init() {
            if (empty()) {
                make_algo();
            }

            std::vector<std::shared_ptr<Transfer>> transfers;
            this->template generate_transfer_operators<IPTruncatedTransfer>(transfers);
            algorithm_->set_transfer_operators(transfers);
            return true;
        }

        bool apply(const Vector &rhs, Vector &x) override {
            assert(algorithm_);
            if (!algorithm_) return false;

            algorithm_->set_box_constraints(this->get_box_constraints());
            return algorithm_->apply(rhs, x);
        }

        inline bool empty() const { return !algorithm_; }

    private:
        std::unique_ptr<MonotoneMultigrid> algorithm_;
        bool is_initialized_{false};
        bool use_patch_smoother_{false};

        void make_algo() {
            InputParameters params;
            params.set("block_size", this->block_size());

            std::shared_ptr<QPSmoother> fine_smoother;
            std::shared_ptr<Smoother> coarse_smoother;

            fine_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            coarse_smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();

            auto direct_solver = std::make_shared<Factorization<Matrix, Vector>>();

            fine_smoother->read(params);
            coarse_smoother->read(params);
            algorithm_ = utopia::make_unique<MonotoneMultigrid>(fine_smoother, coarse_smoother, direct_solver);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_HPP

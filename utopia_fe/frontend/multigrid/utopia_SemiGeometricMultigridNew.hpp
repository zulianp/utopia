#ifndef UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_NEW_HPP
#define UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_NEW_HPP

#include "utopia.hpp"
#include "utopia_Base.hpp"

#include "utopia_fe_Core.hpp"

#include "utopia_SemiGeometricHierarchy.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class FunctionSpace>
    class SemiGeometricMultigridNew
        : public IterativeSolver<typename Traits<FunctionSpace>::Matrix, typename Traits<FunctionSpace>::Vector>,
          public SemiGeometricHierarchy<FunctionSpace> {
    public:
        using Vector = typename Traits<FunctionSpace>::Vector;
        using Matrix = typename Traits<FunctionSpace>::Matrix;
        using Scalar = typename Traits<FunctionSpace>::Scalar;
        using SizeType = typename Traits<FunctionSpace>::SizeType;
        using IndexSet = typename Traits<FunctionSpace>::IndexSet;
        using Comm = typename Traits<FunctionSpace>::Communicator;

        using Super = utopia::IterativeSolver<Matrix, Vector>;
        using IPTransfer = utopia::IPTransfer<Matrix, Vector>;
        using Transfer = utopia::Transfer<Matrix, Vector>;
        using Multigrid = utopia::Multigrid<Matrix, Vector>;

        void read(Input &in) override {
            Super::read(in);

            this->read_hierarchy(in);

            const int n_levels = this->n_coarse_spaces();

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
            UTOPIA_TRACE_REGION_BEGIN("SemiGeometricMultigridNew::update");

            assert(algorithm_);

            if (!is_initialized_) {
                init();
                is_initialized_ = true;
            }

            algorithm_->update(op);

            UTOPIA_TRACE_REGION_END("SemiGeometricMultigridNew::update");
        }

        bool init() {
            if (empty()) {
                make_algo();
            }

            std::vector<std::shared_ptr<Transfer>> transfers;
            this->template generate_transfer_operators<IPTransfer>(transfers);
            algorithm_->set_transfer_operators(transfers);
            return true;
        }

        bool apply(const Vector &rhs, Vector &x) override {
            UTOPIA_TRACE_REGION_BEGIN("SemiGeometricMultigridNew::apply");

            assert(algorithm_);
            bool ok = false;
            if (algorithm_) {
                ok = algorithm_->apply(rhs, x);
            }

            UTOPIA_TRACE_REGION_END("SemiGeometricMultigridNew::apply");
            return ok;
        }

        inline bool empty() const { return !algorithm_; }

    private:
        std::unique_ptr<Multigrid> algorithm_;
        bool is_initialized_{false};

        void make_algo() {
            InputParameters params;
            params.set("block_size", this->block_size());

            auto smoother = std::make_shared<ProjectedGaussSeidel<Matrix, Vector>>();
            auto direct_solver = std::make_shared<KSPSolver<Matrix, Vector>>();

            smoother->read(params);
            algorithm_ = utopia::make_unique<Multigrid>(smoother, direct_solver);
        }
    };

}  // namespace utopia

#endif  // UTOPIA_MONOTONE_SEMI_GEOMETRIC_MULTIGRID_NEW_HPP

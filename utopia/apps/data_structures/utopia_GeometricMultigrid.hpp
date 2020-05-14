#ifndef UTOPIA_GEOMETRIC_MULTIGRID_HPP
#define UTOPIA_GEOMETRIC_MULTIGRID_HPP

#include "utopia_IPTransfer.hpp"
#include "utopia_Input.hpp"
#include "utopia_Multigrid.hpp"
#include "utopia_make_unique.hpp"

#include <memory>

namespace utopia {

    // FIXME complete the overriding process
    template <class FunctionSpace>
    class GeometricMultigrid final : public IterativeSolver<typename FunctionSpace::Matrix,
                                                            typename FunctionSpace::Vector> /*, public Configurable */ {
    public:
        using Matrix = typename FunctionSpace::Matrix;
        using Vector = typename FunctionSpace::Vector;
        using Scalar = typename FunctionSpace::Scalar;
        using SizeType = typename FunctionSpace::SizeType;

        typedef utopia::LinearSolver<Matrix, Vector> Solver;
        // typedef utopia::IterativeSolver<Matrix, Vector>     IterativeSolver;
        typedef utopia::IterativeSolver<Matrix, Vector> Smoother;

        GeometricMultigrid(const std::shared_ptr<Smoother> &smoother,
                           const std::shared_ptr<Solver> &coarse_solver,
                           const bool use_petsc_mg = false) {
            if (use_petsc_mg) {
                algebraic_mg_ =
                    std::make_shared<Multigrid<Matrix, Vector, PETSC_EXPERIMENTAL>>(smoother, coarse_solver);
            } else {
                algebraic_mg_ = std::make_shared<Multigrid<Matrix, Vector>>(smoother, coarse_solver);
            }

            algebraic_mg_->must_generate_masks(false);
        }

        void read(Input &in) override { algebraic_mg_->read(in); }

        void update(const std::shared_ptr<const Matrix> &op) override { algebraic_mg_->update(op); }

        bool apply(const Vector &rhs, Vector &sol) override { return algebraic_mg_->apply(rhs, sol); }

        GeometricMultigrid *clone() const override {
            auto gmg = utopia::make_unique<GeometricMultigrid>(*this);
            return gmg.release();
        }

        inline void verbose(const bool &val) override { algebraic_mg_->verbose(val); }

        inline void max_it(const SizeType &it) override { algebraic_mg_->max_it(it); }

        inline SizeType max_it() const override { return algebraic_mg_->max_it(); }

        inline void atol(const double &tol) override { algebraic_mg_->atol(tol); }

        bool init(FunctionSpace &space, const SizeType n_levels) { return init(make_ref(space), n_levels); }

        bool init(const std::shared_ptr<FunctionSpace> &space, const SizeType n_levels) {
            if (n_levels < 2) {
                std::cerr << "n_levels must be at least 2" << std::endl;
                return false;
            }

            spaces_.resize(n_levels);
            spaces_[0] = space;

            std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> transfers(n_levels - 1);

            for (SizeType i = 1; i < n_levels; ++i) {
                spaces_[i] = spaces_[i - 1]->uniform_refine();

                auto I = std::make_shared<Matrix>();
                spaces_[i - 1]->create_interpolation(*spaces_[i], *I);
                assert(!empty(*I));
                transfers[i - 1] = std::make_shared<IPTransfer<Matrix, Vector>>(I);
            }

            algebraic_mg_->set_transfer_operators(transfers);
            return true;
        }

        bool init(const std::vector<std::shared_ptr<FunctionSpace>> &spaces) {
            const SizeType n_levels = spaces.size();

            if (n_levels < 2) {
                std::cerr << "n_levels must be at least 2" << std::endl;
                return false;
            }

            spaces_ = spaces;

            std::vector<std::shared_ptr<Transfer<Matrix, Vector>>> transfers(n_levels - 1);

            for (SizeType i = 1; i < n_levels; ++i) {
                auto I = std::make_shared<Matrix>();
                spaces_[i - 1]->create_interpolation(*spaces_[i], *I);
                assert(!empty(*I));
                transfers[i - 1] = std::make_shared<IPTransfer<Matrix, Vector>>(I);
            }

            algebraic_mg_->set_transfer_operators(transfers);
            return true;
        }

        FunctionSpace &fine_space() { return *spaces_.back(); }

        const FunctionSpace &fine_space() const { return *spaces_.back(); }

        std::shared_ptr<FunctionSpace> fine_space_ptr() { return spaces_.back(); }

        std::shared_ptr<const FunctionSpace> fine_space_ptr() const { return spaces_.back(); }

    private:
        std::shared_ptr<LinearMultiLevel<Matrix, Vector>> algebraic_mg_;
        std::vector<std::shared_ptr<FunctionSpace>> spaces_;
    };

}  // namespace utopia

#endif  // UTOPIA_GEOMETRIC_MULTIGRID_HPP

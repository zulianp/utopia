#ifndef UTOPIA_SEMISMOOTH_NEWTON_HPP
#define UTOPIA_SEMISMOOTH_NEWTON_HPP

#include "utopia_Core.hpp"
#include "utopia_QPSolver.hpp"

#include "utopia_polymorphic_LinearSolver.hpp"

#include <memory>
#include <vector>

namespace utopia {

    template <class Matrix, class Vector, int Backend = Traits<Vector>::Backend>
    class SemismoothNewton final : public QPSolver<Matrix, Vector> {
        using Scalar = typename Traits<Vector>::Scalar;
        using SizeType = typename Traits<Vector>::SizeType;
        using Layout = typename Traits<Vector>::Layout;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using Super = utopia::QPSolver<Matrix, Vector>;

    public:
        inline void set_linear_solver(const std::shared_ptr<LinearSolver> &linear_solver) {
            linear_solver_ = linear_solver;
        }

        explicit SemismoothNewton(
            std::shared_ptr<LinearSolver> linear_solver = std::make_shared<OmniLinearSolver<Matrix, Vector>>());
        ~SemismoothNewton() override;
        SemismoothNewton *clone() const override;

        void read(Input &in) override;
        void print_usage(std::ostream &os) const override;

        bool smooth(const Vector &, Vector &) override;
        bool apply(const Vector &b, Vector &x) override;
        void init_memory(const Layout &layout) override;
        void update(const std::shared_ptr<const Matrix> &op) override;
        void fallback_solver(const std::shared_ptr<LinearSolver> &fallback_solver);

    private:
        class Details;

        std::shared_ptr<LinearSolver> linear_solver_;
        std::unique_ptr<Details> details_;
    };

}  // namespace utopia

#endif  // UTOPIA_SEMISMOOTH_NEWTON_HPP

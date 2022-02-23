#ifndef UTOPIA_POLYMORPHIC_QP_SOLVER_HPP
#define UTOPIA_POLYMORPHIC_QP_SOLVER_HPP

#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_QPSolver.hpp"

#include <memory>
#include <string>

namespace utopia {
    template <class Matrix, class Vector>
    class QPSolverRegistry;

    template <class Matrix, class Vector>
    class OmniQPSolver : public QPSolver<Matrix, Vector> {
    public:
        using Super = utopia::QPSolver<Matrix, Vector>;
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        using LinearSolver = utopia::LinearSolver<Matrix, Vector>;
        using BoxConstraints = utopia::BoxConstraints<Vector>;
        using QPSolver = utopia::QPSolver<Matrix, Vector>;

        using QPSolverRegistry = utopia::QPSolverRegistry<Matrix, Vector>;

    public:
        OmniQPSolver();
        ~OmniQPSolver() override;
        OmniQPSolver *clone() const override;
        bool apply(const Vector &rhs, Vector &sol) override;
        void read(Input &in) override;

        void set(const std::string &backend, const std::string &type);

        void atol(const Scalar &atol_in) override;
        void rtol(const Scalar &rtol_in) override;
        void stol(const Scalar &stol_in) override;
        void max_it(const SizeType &max_it_in) override;
        void verbose(const bool &verbose_in) override;

        void set_selection(const std::shared_ptr<Vector> &selection) override;

        static QPSolverRegistry &registry();

    private:
        std::unique_ptr<QPSolver> impl_;
    };

}  // namespace utopia

#endif  // UTOPIA_POLYMORPHIC_QP_SOLVER_HPP

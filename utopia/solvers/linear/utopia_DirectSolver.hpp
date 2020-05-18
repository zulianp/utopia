#ifndef UTOPIA_DIRECT_SOLVER_HPP
#define UTOPIA_DIRECT_SOLVER_HPP

#include <algorithm>
#include <iostream>
#include "utopia_LinearSolver.hpp"
#include "utopia_Preconditioner.hpp"

namespace utopia {
    template <typename Matrix, typename Vector>
    class DirectSolver : public LinearSolver<Matrix, Vector> {
    public:
        bool apply(const Vector &rhs, Vector &sol) override {
            std::cerr << "[Warning] Subclass of DirectSolver did not override apply(.,.). Using solve..." << std::endl;
            std::cerr << "[Task]    You need to fix your direct solver and override apply(.,.) instead of solve(.,.,.)."
                      << std::endl;
            assert(false && "override this method");
            return this->solve(*this->get_operator(), rhs, sol);
        }

        void read(Input & /*in*/) override {}

        void print_usage(std::ostream & /*os*/) const override {}
    };

}  // namespace utopia

#endif  // UTOPIA_DIRECT_SOLVER_HPP

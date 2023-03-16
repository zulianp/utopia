#ifndef UTOPIA_AMESOS2_SOLVERS_HPP
#define UTOPIA_AMESOS2_SOLVERS_HPP

#include "utopia_Base.hpp"
#include "utopia_DirectSolver.hpp"

#ifdef UTOPIA_WITH_TRILINOS_AMESOS2
#include "Amesos2_config.h"

namespace utopia {
    /**@ingroup     Linear
     * @brief       Class provides interface to Trilinos Amesos2 solvers \n
     *              For setting up basic parameters, one can use classic Amesos2
     * runtime options
     */
    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class Amesos2Solver {};

    template <typename Matrix, typename Vector>
    class Amesos2Solver<Matrix, Vector, TRILINOS> final : public DirectSolver<Matrix, Vector> {
    public:
        Amesos2Solver(const Amesos2Solver &other);
        Amesos2Solver() = delete;
        Amesos2Solver(const std::string &solver_type);
        ~Amesos2Solver() override;

        void update(const std::shared_ptr<const Matrix> &op) override;
        bool apply(const Vector &rhs, Vector &lhs) override;

        void print_usage(std::ostream &os) const override;
        void read(Input &is) override;

        Amesos2Solver *clone() const override;

        int get_nnzLU() const;
        int get_num_preorder() const;
        int get_num_sym_fact() const;
        int get_num_numeric_fact() const;
        int get_num_solve() const;
        bool get_preordering_done() const;
        bool get_sym_factorization_done() const;
        bool get_num_factorization_done() const;

    private:
        bool keep_symbolic_factorization{false};

        const std::string solver_type_;

        class Impl;
        std::unique_ptr<Impl> impl_;

        bool preordering();
        bool num_factorization();
        bool sym_factorization();
    };

}  // namespace utopia

#endif  // UTOPIA_WITH_TRILINOS_AMESOS2
#endif  // UTOPIA_AMESOS2_SOLVERS_HPP
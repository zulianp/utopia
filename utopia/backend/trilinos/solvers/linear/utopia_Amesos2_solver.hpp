#ifndef UTOPIA_AMESOS2_SOLVERS_HPP
#define UTOPIA_AMESOS2_SOLVERS_HPP

#include "utopia_Base.hpp"
#ifdef WITH_TRILINOS_AMESOS2

#include "Amesos2_config.h"

#ifdef HAVE_AMESOS2_KOKKOS

#include "utopia_DirectSolver.hpp"
#include "utopia_PreconditionedSolver.hpp"
#include "utopia_trilinos_LinearSolverFactory.hpp"

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
        typedef UTOPIA_SCALAR(Vector) Scalar;
        typedef UTOPIA_SIZE_TYPE(Vector) SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;

        Amesos2Solver();
        Amesos2Solver(const Amesos2Solver &other);
        ~Amesos2Solver() override;

        void update(const std::shared_ptr<const Matrix> &op) override;
        bool apply(const Vector &rhs, Vector &lhs) override;

        void read(Input &is) override;
        void print_usage(std::ostream &os) const override;

        int get_nnzLU() const;
        int get_num_preorder() const;
        int get_num_sym_fact() const;
        int get_num_numeric_fact() const;
        int get_num_solve() const;
        bool get_preordering_done() const;
        bool get_sym_factorization_done() const;
        bool get_num_factorization_done() const;

        /**
         * @brief      Reads the xml file based on different layout than read
         *
         * @param[in]  path  location of the xml file
         */
        void read_xml(const std::string &path);

        /**
         * @brief      Sets the parameters.
         */
        void check_parameters();  // override;

        Amesos2Solver *clone() const override;

    private:
        class Impl;
        std::unique_ptr<Impl> impl_;

        // bool set_problem();
        // bool set_problem(Matrix &A);
        bool preordering();
        bool num_factorization();
        bool sym_factorization();
    };

}  // namespace utopia

#endif  // HAVE_AMESOS2_KOKKOS
#endif  // UTOPIA_AMESOS2_SOLVERS_HPP
#endif  // WITH_TRILINOS_AMESOS2
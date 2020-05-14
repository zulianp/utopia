#ifndef UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP
#define UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP

#include "utopia_QPSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_LinearSolver.hpp"

#include <memory>

namespace utopia {

    template<class Matrix, class Vector>
    class TaoQPSolver final : public QPSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;
        typedef utopia::LinearSolver<Matrix, Vector> LinearSolver;
        typedef utopia::QPSolver<Matrix, Vector> 	 Super;
        using BoxConstraints = utopia::BoxConstraints<Vector>;

    public:

        TaoQPSolver(const std::shared_ptr<LinearSolver> &linear_solver = std::make_shared<Factorization<Matrix, Vector>>());
        ~TaoQPSolver() override;
        TaoQPSolver * clone() const override;
        bool apply(const Vector &rhs, Vector &sol) override;
        void set_linear_solver(const std::shared_ptr<LinearSolver> &linear_solver);
        void tao_type(const std::string &type);
        // void set_ksp_types(const std::string &ksp_type, const std::string &pc_type, const std::string &solver_package);
        void read(Input &in) override;


        Scalar atol() const override;
        Scalar rtol() const override;
        Scalar stol() const override;

        SizeType max_it() const override;
        bool verbose() const override;
   
        void atol(const Scalar &atol) override;
        void rtol(const Scalar &rtol) override;
        void stol(const Scalar &stol) override;
        void max_it(const SizeType & max_it) override;
        void verbose(const bool &verbose) override;
    private:
        class Impl;
        std::unique_ptr<Impl> impl_;

    };

}

#endif //UTOPIA_PETSC_TAO_QP_SOLVEWR_HPP

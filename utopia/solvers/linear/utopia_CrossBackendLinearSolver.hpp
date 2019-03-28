#ifndef CROSS_BACKEND_LINEAR_SOLVER_HPP
#define CROSS_BACKEND_LINEAR_SOLVER_HPP

#include "utopia_LinearSolver.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Smoother.hpp"
#include "utopia_Conversions.hpp"

namespace utopia {

    template<
        class Matrix,
        class Vector,
        class WantedMatrix,
        class WantedVector,
        class Solver
    >
    class CrossBackendLinearSolver : public LinearSolver<Matrix, Vector> {
    public:
        virtual ~CrossBackendLinearSolver() {}

        virtual bool apply(const Vector &rhs, Vector &sol) override
        {
            backend_convert(rhs, rhs_buff_);
            backend_convert(sol, sol_buff_);

            bool ok = solver_.apply(rhs_buff_, sol_buff_);

            backend_convert(sol_buff_, sol);
            return ok;
        }

        virtual void update(const std::shared_ptr<const Matrix> &op) override {
            if(!mat_buff_) {
                mat_buff_ = std::make_shared<WantedMatrix>();
            }

            backend_convert_sparse(*op, *mat_buff_);
            solver_.update(mat_buff_);
        }

        CrossBackendLinearSolver * clone() const override
        {
            return new CrossBackendLinearSolver();
        }

    private:
        Solver solver_;
        std::shared_ptr<WantedMatrix> mat_buff_;
        WantedVector rhs_buff_;
        WantedVector sol_buff_;
    };


    template<
        class Matrix,
        class Vector,
        class WantedMatrix,
        class WantedVector,
        class Solver
    >
    class CrossBackendLinearSolverAndSmoother : public IterativeSolver<Matrix, Vector>, public Smoother<Matrix, Vector> {
    public:
        virtual ~CrossBackendLinearSolverAndSmoother() {}

        virtual bool apply(const Vector &rhs, Vector &sol) override
        {
            backend_convert(rhs, rhs_buff_);
            backend_convert(sol, sol_buff_);

            bool ok = solver_.apply(rhs_buff_, sol_buff_);

            backend_convert(sol_buff_, sol);
            return ok;
        }

        virtual bool smooth(const Vector &rhs, Vector &x) override
        {
            backend_convert(rhs, rhs_buff_);
            backend_convert(x, sol_buff_);

            bool ok = solver_.smooth(rhs_buff_, sol_buff_);

            backend_convert(sol_buff_, x);
            return ok;
        }

        virtual void update(const std::shared_ptr<const Matrix> &op) override {
            if(!mat_buff_) {
                mat_buff_ = std::make_shared<WantedMatrix>();
            }

            backend_convert_sparse(*op, *mat_buff_);
            solver_.update(mat_buff_);
        }

        CrossBackendLinearSolverAndSmoother * clone() const override
        {
            return new CrossBackendLinearSolverAndSmoother();
        }

        void read(Input &in) override {
             solver_.read(in);
        }

        void print_usage(std::ostream &os = std::cout) const override {
              solver_.print_usage(os);
        }

    private:
        Solver solver_;
        std::shared_ptr<WantedMatrix> mat_buff_;
        WantedVector rhs_buff_;
        WantedVector sol_buff_;
    };
}

#endif //CROSS_BACKEND_LINEAR_SOLVER_HPP

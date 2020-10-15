#ifndef UTOPIA_SP_BLOCK_CONJUGATE_GRADIENT_HPP
#define UTOPIA_SP_BLOCK_CONJUGATE_GRADIENT_HPP

#include "utopia_ForwardDeclarations.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia_Layout.hpp"
#include "utopia_LinearSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_Traits.hpp"

#include <iostream>

namespace utopia {

    /**
    @brief Solves the following saddle-point problem:
        | A_m		0		B_t | | sol_m | = | rhs_m |
        | 0			A_s		D_t | | sol_s | = | rhs_s |
        | B 		D 		0   | | lagr  | = | 0     |
    m := master, s := slave
    It uses the schur complement
    Warning: rows with boundary conditions in B_t and D_t have to be handled outside.

    */
    template <class Matrix, class Vector, int Backend = Traits<Matrix>::Backend>
    class SPBlockConjugateGradient {
    public:
        using Scalar = UTOPIA_SCALAR(Vector);

        SPBlockConjugateGradient()
            : op_m(std::make_shared<Factorization<Matrix, Vector>>()),
              op_s(std::make_shared<Factorization<Matrix, Vector>>()),
              atol_(1e-8),
              rtol_(1e-8) {}

        void set_master_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &op_m) { this->op_m = op_m; }

        void set_slave_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &op_s) { this->op_s = op_s; }

        void set_prec_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &solver) {
            if (prec_) {
                prec_->set_linear_solver(solver);
            } else {
                std::cerr << "[Error] set preconditioner before calling set_prec_solver" << std::endl;
                assert(prec_);
            }
        }

        void set_master_sweeps(const int n_sweeps) { master_sweeps_ = n_sweeps; }

        void set_master_max_it(const int n_its) { master_max_it_ = n_its; }

        inline void verbose(const bool verbose) { verbose_ = verbose; }

        inline void atol(const Scalar atol) { atol_ = atol; }

        inline void max_it(const int max_it) { max_it_ = max_it; }

        // inline void rtol(const Scalar rtol)
        // {
        // 	rtol_ = rtol;
        // }

        void update(const std::shared_ptr<Matrix> &A_m,
                    const std::shared_ptr<Matrix> &A_s,
                    const std::shared_ptr<Matrix> &B,
                    const std::shared_ptr<Matrix> &D,
                    const std::shared_ptr<Matrix> &B_t,
                    const std::shared_ptr<Matrix> &D_t) {
            set_up(B, D, B_t, D_t);
            update(A_m, A_s);
        }

        void set_up(const std::shared_ptr<Matrix> &B,
                    const std::shared_ptr<Matrix> &D,
                    const std::shared_ptr<Matrix> &B_t,
                    const std::shared_ptr<Matrix> &D_t) {
            this->B = B;
            this->B_t = B_t;
            this->D = D;
            this->D_t = D_t;
        }

        void update(const std::shared_ptr<Matrix> &A_m, const std::shared_ptr<Matrix> &A_s) {
            this->A_m = A_m;
            this->A_s = A_s;

            op_m->update(A_m);
            op_s->update(A_s);

            set_iterations_master(master_sweeps_);

            if (prec_) {
                prec_->update(*this);
            }
        }

        bool apply(const Vector &rhs_m, const Vector &rhs_s, Vector &sol_m, Vector &sol_s, Vector &lagr) {
            if (prec_) {
                return preconditioned_solve(rhs_m, rhs_s, sol_m, sol_s, lagr);
            } else {
                return unpreconditioned_solve(rhs_m, rhs_s, sol_m, sol_s, lagr);
            }
        }

        // http://ta.twi.tudelft.nl/nw/users/vuik/talks/norwich_2014.pdf
        // void init_simple_preconditioner()
        // {
        // 	Matrix d_m = diag(1./diag(*A_m));
        // 	Matrix d_s = diag(1./diag(*A_s));

        // 	P = std::make_shared<Matrix>();
        // 	*P = (*B) * d_m * transpose(*B) + (*D) * d_s * transpose(*D);

        // 	if(!op_prec) {
        // 		op_prec = std::make_shared<Factorization<Matrix, Vector>>();
        // 	}

        // 	op_prec->update(P);
        // }

        class BlockPreconditioner {
        public:
            virtual ~BlockPreconditioner() = default;
            virtual bool apply(const Vector &r, Vector &z) const = 0;
            virtual void update(SPBlockConjugateGradient &solver) = 0;
            virtual void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &solver) = 0;
        };

        // http://ta.twi.tudelft.nl/nw/users/vuik/talks/norwich_2014.pdf
        class SIMPLEPreconditioner final : public BlockPreconditioner {
        public:
            bool apply(const Vector &r, Vector &z) const override { return op_prec->apply(r, z); }

            void update(SPBlockConjugateGradient &solver) override {
                auto &D = *solver.D;
                auto &B = *solver.B;

                auto &D_t = *solver.D_t;
                auto &B_t = *solver.B_t;

                auto &A_m = *solver.A_m;
                auto &A_s = *solver.A_s;

                Matrix d_m = diag(1. / diag(A_m));
                Matrix d_s = diag(1. / diag(A_s));

                P = std::make_shared<Matrix>();
                *P = B * d_m * B_t + D * d_s * D_t;

                if (!op_prec) {
                    op_prec = std::make_shared<Factorization<Matrix, Vector>>();
                }

                op_prec->update(P);
            }

            inline void set_linear_solver(const std::shared_ptr<LinearSolver<Matrix, Vector>> &solver) override {
                op_prec = solver;
            }

        private:
            std::shared_ptr<Matrix> P;
            std::shared_ptr<LinearSolver<Matrix, Vector>> op_prec;
        };

        // class KlawonWildlundPreconditioner final : public BlockPreconditioner {
        // public:
        // 	bool apply(const Vector &r, Vector &z) const override
        // 	{
        // 		// op_prec->apply(r, z);
        // 	}

        // 	void update(SPBlockConjugateGradient &solver) override
        // 	{
        // 		auto &D = *solver.D;
        // 		auto &B = *solver.B;
        // 		auto &A_m = *solver.A_m;
        // 		auto &A_s = *solver.A_s;

        // 		Matrix d_m = diag(1./diag(A_m));
        // 		Matrix d_s = diag(1./diag(A_s));

        // 		P = std::make_shared<Matrix>();
        // 		*P = (*B) * transpose(*B) + (*D) * transpose(*D);

        // 		if(!op_prec) {
        // 			op = std::make_shared<Factorization<Matrix, Vector>>();
        // 		}

        // 		op->update(P);
        // 	}

        // private:
        // 	std::shared_ptr<Matrix> P, BDtxBD, BDxBDt;
        // 	std::shared_ptr< LinearSolver<Matrix, Vector> > op;
        // };

        // void use_klawon_wildlund_preconditioner()
        // {

        // }

        void use_simple_preconditioner() { prec_ = std::make_shared<SIMPLEPreconditioner>(); }

    private:
        std::shared_ptr<Matrix> A_m, A_s;

        std::shared_ptr<Matrix> D, D_t;
        std::shared_ptr<Matrix> B, B_t;

        std::shared_ptr<Matrix> BDxBDt;
        std::shared_ptr<Matrix> BDtxBD;

        std::shared_ptr<LinearSolver<Matrix, Vector>> op_m;
        std::shared_ptr<LinearSolver<Matrix, Vector>> op_s;

        Vector buff_m, solved_m, buff_s, solved_s, r;

        bool verbose_{false};
        Scalar atol_, rtol_;
        int max_it_{1000};

        int master_sweeps_{-1};
        int master_max_it_{-1};

        std::shared_ptr<BlockPreconditioner> prec_;

        Vector p, q, Ap, r_new, z, z_new;

        void apply_preconditioner(const Vector &r, Vector &z) {
            assert(prec_);
            prec_->apply(r, z);
        }

        inline void apply_op(const Vector &p, Vector &Ap) {
            buff_m = (*B_t) * p;
            buff_s = (*D_t) * p;

            if (empty(solved_m) || size(buff_m) != size(solved_m)) {
                solved_m.zeros(layout(buff_m));
            } else {
                solved_m.set(0.);
            }

            if (empty(solved_s) || size(buff_s) != size(solved_s)) {
                solved_s.zeros(layout(buff_s));
            } else {
                solved_s.set(0.);
            }

            op_m->apply(buff_m, solved_m);
            op_s->apply(buff_s, solved_s);

            Ap = (*B) * solved_m + (*D) * solved_s;
        }

        inline void residual(const Vector &rhs_m, const Vector &rhs_s, const Vector &p, Vector &r) {
            buff_m = rhs_m - (*B_t) * p;
            buff_s = rhs_s - (*D_t) * p;

            if (empty(solved_m) || size(buff_m) != size(solved_m)) {
                solved_m.zeros(layout(buff_m));
            } else {
                solved_m.set(0.);
            }

            if (empty(solved_s) || size(buff_s) != size(solved_s)) {
                solved_s.zeros(layout(buff_s));
            } else {
                solved_s.set(0.);
            }

            op_m->apply(buff_m, solved_m);
            op_s->apply(buff_s, solved_s);

            r = (*B) * solved_m + (*D) * solved_s;
        }

        bool preconditioned_solve(const Vector &rhs_m,
                                  const Vector &rhs_s,
                                  Vector &sol_m,
                                  Vector &sol_s,
                                  Vector &lagr) {
            if (empty(lagr) || size(lagr).get(0) != size(*B).get(0)) {
                lagr.zeros(row_layout(*B));
            }

            this->residual(rhs_m, rhs_s, lagr, r);

            Scalar it = 0;
            Scalar beta = 0., alpha = 1., r_norm = 9e9;

            Vector lagr_old = lagr;

            z.zeros(layout(r));
            z_new.zeros(layout(r));

            this->apply_preconditioner(r, z);

            p = z;

            bool converged = false;
            while (!converged) {
                this->apply_op(p, Ap);

                alpha = dot(r, z) / dot(p, Ap);
                lagr += alpha * p;
                r_new = r - alpha * Ap;

                r_norm = norm2(r_new);

                if (verbose_) {
                    Scalar diff_norm = norm2(lagr - lagr_old);
                    utopia::out() << it << " " << r_norm << " " << diff_norm << std::endl;
                }

                if (r_norm <= atol_) {
                    converged = true;
                    break;
                }

                this->apply_preconditioner(r_new, z_new);

                beta = dot(z_new, r_new) / dot(z, r);

                p = z_new + beta * p;
                r = r_new;
                z = z_new;

                it++;

                if (!converged && it > max_it_) {
                    break;
                }

                lagr_old = lagr;
            }

            if (empty(sol_m) || size(sol_m) != size(rhs_m)) {
                sol_m.zeros(layout(rhs_m));
            }

            if (empty(sol_s) || size(sol_s) != size(rhs_s)) {
                sol_s.zeros(layout(rhs_s));
            }

            op_m->apply(rhs_m - (*B_t) * lagr, sol_m);
            op_s->apply(rhs_s - (*D_t) * lagr, sol_s);

            return converged;
        }

        bool unpreconditioned_solve(const Vector &rhs_m,
                                    const Vector &rhs_s,
                                    Vector &sol_m,
                                    Vector &sol_s,
                                    Vector &lagr) {
            if (empty(lagr) || size(lagr).get(0) != size(*B).get(0)) {
                lagr.zeros(row_layout(*B));
            }

            this->residual(rhs_m, rhs_s, lagr, r);

            bool converged = false;

            Scalar it = 0;
            Scalar rho = 1., rho_1 = 1., beta = 0., alpha = 1., r_norm = 9e9;

            Vector lagr_old = lagr;

            while (!converged) {
                rho = dot(r, r);

                if (rho == 0.) {
                    converged = true;
                    break;
                }

                if (it > 0) {
                    beta = rho / rho_1;
                    p = r + beta * p;
                } else {
                    p = r;
                }

                this->apply_op(p, q);
                alpha = rho / dot(p, q);

                lagr += alpha * p;
                r -= alpha * q;

                rho_1 = rho;
                it++;
                r_norm = norm2(lagr - lagr_old);

                converged = r_norm < atol_;
                it++;

                if (!converged && it >= max_it_) {
                    break;
                }

                if (verbose_) {
                    utopia::out() << it << " " << r_norm << std::endl;
                }

                lagr_old = lagr;
            }

            if (empty(sol_m) || size(sol_m) != size(rhs_m)) {
                sol_m.zeros(layout(rhs_m));
            }

            if (empty(sol_s) || size(sol_s) != size(rhs_s)) {
                sol_s.zeros(layout(rhs_s));
            }

            set_iterations_master(master_max_it_);

            op_m->apply(rhs_m - (*B_t) * lagr, sol_m);
            op_s->apply(rhs_s - (*D_t) * lagr, sol_s);

            return converged;
        }

        void set_iterations_master(const int iterations) {
            if (iterations > 0) {
                auto op_m_iterative = std::dynamic_pointer_cast<IterativeSolver<Matrix, Vector>>(op_m);
                if (op_m_iterative) {
                    op_m_iterative->max_it(iterations);
                }
            }
        }
    };

}  // namespace utopia

#endif  // UTOPIA_SP_BLOCK_CONJUGATE_GRADIENT_HPP

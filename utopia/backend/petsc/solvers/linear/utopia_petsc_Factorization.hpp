#ifndef UTOPIA_FACTORIZATION_HPP
#define UTOPIA_FACTORIZATION_HPP

#include "utopia_DirectSolver.hpp"
#include "utopia_LinearSolverInterfaces.hpp"
#include "utopia_SolverType.hpp"
#include "utopia_petsc_KSPSolver.hpp"

namespace utopia {

    //////////////////////////////////////////////////////////
    template <typename Matrix, typename Vector>
    class Factorization<Matrix, Vector, PETSC> : public DirectSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

    public:
        Factorization() {
            strategy_.ksp_type(KSPPREONLY);
            strategy_.set_initial_guess_non_zero(false);
            strategy_.pc_type(PCLU);
            strategy_.solver_package(default_package());
        }

        inline constexpr static SolverPackage default_package() {
#ifdef PETSC_HAVE_SUPERLU_DIST
            return MATSOLVERSUPERLU_DIST;
#else  // PETSC_HAVE_SUPERLU_DIST
#ifdef PETSC_HAVE_MUMPS
            return MATSOLVERMUMPS;
#else  // PETSC_HAVE_MUMPS

#ifdef PETSC_HAVE_SUPERLU
            return MATSOLVERSUPERLU;
#else   // PETSC_HAVE_SUPERLU
            return MATSOLVERPETSC;
#endif  // PETSC_HAVE_SUPERLU
#endif  // PETSC_HAVE_MUMPS
#endif  // PETSC_HAVE_SUPERLU_DIST
        }

        Factorization(const std::string &sp, const std::string &pct) {
            strategy_.ksp_type(KSPPREONLY);
            strategy_.set_initial_guess_non_zero(false);
            strategy_.pc_type(pct);
            strategy_.solver_package(sp);
        }

        void set_type(const std::string &lib, const std::string &type) {
            strategy_.solver_package(lib);
            strategy_.pc_type(type);
        }

        inline bool apply(const Vector &b, Vector &x) override {
            if (this->is_sparse()) {
                return strategy_.apply(b, x);
            } else {
                return dense_strategy_.apply(b, x);
            }
        }

        inline void update(const std::shared_ptr<const Matrix> &op) override {
            if (op->is_sparse()) {
                is_sparse_ = true;
                strategy_.update(op);
            } else {
                is_sparse_ = false;
                dense_strategy_.update(op);
            }
        }

        Factorization *clone() const override { return new Factorization(*this); }

        void describe(std::ostream &os) const { strategy_.describe(os); }

        inline void read(Input &is) override {
            DirectSolver<Matrix, Vector>::read(is);

            // FIXME check before if it is a factorization
            strategy_.read(is);
        }

    private:
        class Strategy : public KSPSolver<Matrix, Vector> {
        public:
            using KSPSolver<Matrix, Vector>::KSPSolver;
            /*@todo use : PCFactorSetReuseFill(PC pc,PetscBool flag)  // to keep factorization from prev. levels
             */

            void set_ksp_options(KSP &ksp) override {
                this->reset_preconditioner();
                this->ksp_type(KSPPREONLY);
                this->set_initial_guess_non_zero(false);

                KSPSolver<Matrix, Vector>::set_ksp_options(ksp);
            }
        };

        class DenseStrategy {
        public:
            inline void update(const std::shared_ptr<const Matrix> &op) {
                destroy();
                op_ = std::make_shared<Matrix>();
                op_->copy(*op);

                auto &raw_mat = op_->raw_type();

                PetscErrorCode ierr = 0;

                ierr = MatGetOrdering(raw_mat, MATORDERINGNATURAL, &isr, &isc);
                assert(ierr == 0);

                ierr = MatLUFactor(raw_mat, isr, isc, &info);
                assert(ierr == 0);

                empty_ = ierr == 0;
            }

            inline bool apply(const Vector &b, Vector &x) {
                if (empty_) return false;

                auto raw_mat = op_->raw_type();

                if (x.empty()) {
                    x.zeros(layout(b));
                }

                PetscErrorCode ierr = MatSolve(raw_mat, b.raw_type(), x.raw_type());
                return ierr == 0;
            }

            void destroy() {
                if (isr != nullptr) {
                    ISDestroy(&isr);
                    isr = nullptr;
                }
                if (isc != nullptr) {
                    ISDestroy(&isc);
                    isc = nullptr;
                }

                op_ = nullptr;
                empty_ = true;
            }

            ~DenseStrategy() { destroy(); }

        private:
            std::shared_ptr<Matrix> op_;
            IS isr{nullptr}, isc{nullptr};
            MatFactorInfo info;
            bool empty_{true};
        };

        Strategy strategy_;
        DenseStrategy dense_strategy_;
        bool is_sparse_{true};

    public:
        inline Strategy &strategy() { return strategy_; }

        inline bool is_proxy() const { return false; }
        inline bool is_sparse() const { return is_sparse_; }
    };  // namespace utopia
}  // namespace utopia

#endif  // UTOPIA_FACTORIZATION_HPP

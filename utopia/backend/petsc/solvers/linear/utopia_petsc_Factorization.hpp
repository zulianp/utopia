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

        inline bool apply(const Vector &b, Vector &x) override { return strategy_.apply(b, x); }

        inline void update(const std::shared_ptr<const Matrix> &op) override { strategy_.update(op); }

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

        Strategy strategy_;

    public:
        inline Strategy &strategy() { return strategy_; }
    };
}  // namespace utopia

#endif  // UTOPIA_FACTORIZATION_HPP

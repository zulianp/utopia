#ifndef UTOPIA_KSP_MF
#define UTOPIA_KSP_MF
#include <cstring>
#include <string>
#include "utopia_MatrixFreeLinearSolver.hpp"
#include "utopia_petsc_KSPSolver.hpp"

namespace utopia {

    PetscErrorCode UtopiaOperatorApplyShell(Mat m, Vec X, Vec Y);

    template <typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class KSP_MF {};

    // this class could be merged with KSP class - not sure, if we care or not...
    template <typename Matrix, typename Vector>
    class KSP_MF<Matrix, Vector, PETSC> : public OperatorBasedLinearSolver<Matrix, Vector> {
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        typedef utopia::OperatorBasedLinearSolver<Matrix, Vector> OperatorBasedLinearSolver;

        using Preconditioner = utopia::Preconditioner<Vector>;
        typedef utopia::IterativeSolver<Matrix, Vector> IterativeSolver;

        static_assert(Traits<Matrix>::Backend == utopia::PETSC, "utopia::KSP_MF:: only works with petsc types");

    public:
        KSP_MF() : OperatorBasedLinearSolver() {
            ksp_.pc_type("bjacobi");
            ksp_.verbose(true);
        }
        ~KSP_MF() override = default;

        void read(Input &in) override {
            OperatorBasedLinearSolver::read(in);
            ksp_.read(in);
        }

        bool solve(const Operator<Vector> &A, const Vector &b, Vector &x) {
            update(A);
            return apply(b, x);
        }

        void update(const Operator<Vector> &A) {
            ksp_.handle_reset(A.comm());
            ksp_.set_monitor_options(ksp_.implementation());

            Mat op_mat;
            auto A_nonconst_ref = const_cast<Operator<Vector> *>(&A);

            MatCreateShell(A.comm().get(),
                           A.local_size().get(0),
                           A.local_size().get(1),
                           PETSC_DECIDE,
                           PETSC_DECIDE,
                           A_nonconst_ref,
                           &op_mat);

            // TODO:: check if preconditioners are updated correctly
            // in theory they both coud be just op,  or 1 op and 1 matrix

            MatSetType(op_mat, MATSHELL);
            MatShellSetOperation(op_mat, MATOP_MULT, (void (*)(void))UtopiaOperatorApplyShell);
            KSPSetOperators(ksp_.implementation(), op_mat, op_mat);
        }

        bool apply(const Vector &b, Vector &x) override {
            return ksp_.apply(b, x);
            ksp_.solution_status(this->solution_status_);
        }

        void print_usage(std::ostream &os) const override {
            OperatorBasedLinearSolver::print_usage(os);
            ksp_.print_usage();
        }

        KSP_MF<Matrix, Vector, PETSC> *clone() const override { return new KSP_MF<Matrix, Vector, PETSC>(*this); }

        void update(const std::shared_ptr<const Matrix> &op) override { ksp_.update(op); }

        void set_preconditioner(const std::shared_ptr<Preconditioner> &precond) { ksp_.set_preconditioner(precond); }

        // necessary, as ksp_ is resetting options with every apply...
        void atol(const Scalar &atol_in) override { ksp_.atol(atol_in); }
        void stol(const Scalar &stol_in) override { ksp_.stol(stol_in); }
        void rtol(const Scalar &rtol_in) override { ksp_.rtol(rtol_in); }
        void max_it(const SizeType &max_it_in) override { ksp_.max_it(max_it_in); }
        void verbose(const bool &verbose_in) override { ksp_.verbose(verbose_in); }

        void pc_type(const std::string &pc_type_flg) { ksp_.pc_type(pc_type_flg); }
        void ksp_type(const std::string &ksp_type_flg) { ksp_.ksp_type(ksp_type_flg); }
        void solver_package(const std::string &package_flg) { ksp_.solver_package(package_flg); }
        bool smooth(const Vector &rhs, Vector &x) { return ksp_.smooth(rhs, x); }
        void set_monitor_options(KSP &ksp) const { ksp_.set_monitor_options(ksp); }

    private:
        KSPSolver<Matrix, Vector, PETSC> ksp_;
    };

}  // namespace utopia

#endif  // UTOPIA_KSP_MF

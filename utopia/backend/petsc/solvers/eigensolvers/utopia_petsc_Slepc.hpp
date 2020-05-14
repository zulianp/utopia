#ifndef UTOPIA_PETSC_SLEPC_H
#define UTOPIA_PETSC_SLEPC_H

#include "utopia_Core.hpp"
#include "utopia_EigenSolver.hpp"
#include "utopia_Wrapper.hpp"
#include "utopia_Clonable.hpp"
#include "utopia_Traits.hpp"

#include <vector>
#include <string>

#include <slepceps.h>

namespace utopia {

    template<typename Matrix, typename Vector, int Backend = Traits<Matrix>::Backend>
    class SlepcSolver;

    template<typename Matrix, typename Vector>
    class SlepcSolver<Matrix, Vector, PETSC_EXPERIMENTAL> final : public EigenSolver<Matrix, Vector> {
    public:
        using Scalar = typename utopia::Traits<Vector>::Scalar;
        using SizeType = typename utopia::Traits<Vector>::SizeType;

        SlepcSolver(const std::vector<std::string> problem_types = {"hermitian", "non_hermitian", "generalized_hermitian", "generalized_non_hermitian", "generalized_hermitian_SPD_B", "generalized_hermitian_indefinite"},
                    const std::vector<std::string> portions_of_spectrum = {"largest_magnitude", "smallest_magnitude", "largest_real", "smallest_real", "largest_imaginary", "smallest_imaginary", "closest_to_target", "closest_to_target_real", "closest_to_target_imaginary", "all_in_region"},
                    const std::vector<std::string> solver_types = {"krylovschur", "power", "subspace", "arnoldi", "lanczos", "gd", "jd", "rqcg", "lobpcg", "ciss", "lapack", "arpack", "blzpack", "trlan", "blopex", "primme", "feast"});

        ~SlepcSolver() override;

        SlepcSolver * clone() const override;

        void problem_type(const std::string & type);

        const std::string & problem_type() const;

        void solver_type(const std::string & type);

        const std::string & solver_type() const;

        void portion_of_spectrum(const std::string & type) override;

        bool solve(const Matrix & A) override;

        bool solve(const Matrix & A, const Matrix & B) override;

        bool print_eigenpairs() override;

        void get_eigenpairs(const SizeType & i, Scalar & iegr, Scalar & eigi, Vector & vr, Vector & vi) override;
        void get_real_eigenpair(const SizeType & i, Scalar & iegr, Vector & vr) override;

    private:
        void initialize(const MPI_Comm & comm);
        void reinitialize(const MPI_Comm & comm);

    private:
        EPS eps_;
        bool initialized_;
        bool solved_;

        const std::vector<std::string> problem_types_;
        const std::vector<std::string> portions_of_spectrum_;
        const std::vector<std::string> solver_types_;

        std::string portion_of_spectrum_;
        std::string eps_type_;
        std::string problem_type_;
        std::string solver_type_;
    };

}

#endif //UTOPIA_PETSC_SLEPC_H
#include "utopia_FactoryMethod.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_LinearSolverFactory.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"
#include "utopia_polymorphic_QPSolver.hpp"

// includes
// HOMEMADE
#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGaussSeidel.hpp"
#include "utopia_ProjectedGradient.hpp"
#include "utopia_SemismoothNewton.hpp"

// PETSC
#ifdef UTOPIA_WITH_PETSC
#include "petsctao.h"
#include "utopia_petsc_Factorizations.hpp"
#include "utopia_petsc_LinearSolvers.hpp"
#include "utopia_petsc_RowView.hpp"
#include "utopia_petsc_SemismoothNewton.hpp"
#include "utopia_petsc_TaoQPSolver.hpp"

//
#ifdef UTOPIA_WITH_BLAS
#include "utopia_blas.hpp"
#include "utopia_blas_Array.hpp"
#include "utopia_petsc_RASPatchSmoother.hpp"
#endif  // UTOPIA_WITH_BLAS

#endif  // UTOPIA_WITH_PETSC

#include <functional>
#include <map>
#include <memory>
#include <string>

#include "utopia_FactoryMethod.hpp"

namespace utopia {

    template <class Matrix, class Vector>
    class QPSolverRegistry {
    public:
        using Scalar = typename Traits<Vector>::Scalar;
        using QPSolver = utopia::QPSolver<Matrix, Vector>;
        using QPSolverPtr = std::unique_ptr<QPSolver>;
        using LinearSolverPtr = std::unique_ptr<utopia::LinearSolver<Matrix, Vector>>;

        using FactoryMethod_t = utopia::IFactoryMethod<QPSolver>;

        template <class Alg>
        using QPSFactoryMethod = FactoryMethod<QPSolver, Alg>;

        std::map<std::string, std::map<std::string, std::unique_ptr<FactoryMethod_t>>> solvers;

        inline static QPSolverRegistry &instance() {
            static QPSolverRegistry instance_;
            return instance_;
        }

        template <class SolverT>
        inline void register_solver(const std::string &backend, const std::string &type) {
            solvers[backend][type] = utopia::make_unique<QPSFactoryMethod<SolverT>>();
        }

        template <class SolverT>
        inline void register_solver() {
            solvers[SolverT::backend()][SolverT::solver_type()] = utopia::make_unique<QPSFactoryMethod<SolverT>>();
        }

        QPSolverPtr find(const std::string &backend, const std::string &type) const {
            auto b_it = solvers.find(backend);
            if (b_it == solvers.end()) {
                utopia::err() << "[Error] backend " << backend << " not found using default" << std::endl;
                return default_solver();
            }

            auto s_it = b_it->second.find(type);

            if (s_it == b_it->second.end()) {
                utopia::err() << "[Error] solver " << type << " not found using default" << std::endl;
                return default_solver();
            }

            return s_it->second->make();
        }

        inline static LinearSolverPtr default_linear_solver() {
            auto ls = utopia::make_unique<OmniLinearSolver<Matrix, Vector>>();
            InputParameters in;
            in.set("type", Solver::direct());
            ls->read(in);
            return ls;
        }

        inline static QPSolverPtr default_solver() {
            return utopia::make_unique<SemismoothNewton<Matrix, Vector>>(default_linear_solver());
        }

        QPSolverRegistry() {
            {
                // Homemade
                using HomeMadeSemismoothNewton = utopia::SemismoothNewton<Matrix, Vector>;
                using HomeMadeProjectedGradient = utopia::ProjectedGradient<Matrix, Vector>;
                using HomeMadeProjectedConjugateGradient = utopia::ProjectedConjugateGradient<Matrix, Vector>;
                using HomeMadeProjectedGaussSeidel = utopia::ProjectedGaussSeidel<Matrix, Vector>;

                register_solver<HomeMadeSemismoothNewton>("any", "ssnewton");
                register_solver<HomeMadeProjectedGradient>("any", "pg");
                register_solver<HomeMadeProjectedConjugateGradient>("any", "pcg");
                register_solver<HomeMadeProjectedGaussSeidel>("any", "pgs");
            }

#ifdef UTOPIA_WITH_PETSC
            {
                // Petsc
                using PetscSemiSmoothNewton = utopia::SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL>;
                using PetscTaoQPSolver = utopia::TaoQPSolver<Matrix, Vector>;

                register_solver<PetscSemiSmoothNewton>("petsc", "ssnewton");
                register_solver<PetscTaoQPSolver>("petsc", "taoqp");

#ifdef UTOPIA_WITH_BLAS
                using PetscRASPatchSmoother = utopia::RASPatchSmoother<Matrix, utopia::BlasMatrix<Scalar>>;
                register_solver<PetscRASPatchSmoother>("petsc", "patch_smoother");
#endif  // UTOPIA_WITH_BLAS
            }
#endif  // UTOPIA_WITH_PETSC
        }
    };

    template <class Matrix, class Vector>
    void OmniQPSolver<Matrix, Vector>::set(const std::string &backend, const std::string &type) {
        impl_ = QPSolverRegistry::instance().find(backend, type);
    }

    template <class Matrix, class Vector>
    void OmniQPSolver<Matrix, Vector>::read(Input &in) {
        Super::read(in);

        std::string backend = Traits<Vector>::backend_info().get_name();
        std::string type = "ssnewton";

        in.get("backend", backend);
        in.get("type", type);

        impl_ = QPSolverRegistry::instance().find(backend, type);
        impl_->read(in);
    }

    template <class Matrix, class Vector>
    void OmniQPSolver<Matrix, Vector>::atol(const Scalar &atol_in) {
        Super::atol(atol_in);
        if (impl_) {
            impl_->atol(atol_in);
        }
    }

    template <class Matrix, class Vector>
    void OmniQPSolver<Matrix, Vector>::rtol(const Scalar &rtol_in) {
        Super::rtol(rtol_in);
        if (impl_) {
            impl_->rtol(rtol_in);
        }
    }

    template <class Matrix, class Vector>
    void OmniQPSolver<Matrix, Vector>::stol(const Scalar &stol_in) {
        Super::stol(stol_in);
        if (impl_) {
            impl_->stol(stol_in);
        }
    }

    template <class Matrix, class Vector>
    void OmniQPSolver<Matrix, Vector>::max_it(const SizeType &max_it_in) {
        Super::max_it(max_it_in);
        if (impl_) {
            impl_->max_it(max_it_in);
        }
    }

    template <class Matrix, class Vector>
    void OmniQPSolver<Matrix, Vector>::verbose(const bool &verbose_in) {
        Super::verbose(verbose_in);
        if (impl_) {
            impl_->verbose(verbose_in);
        }
    }

    template <class Matrix, class Vector>
    QPSolverRegistry<Matrix, Vector> &OmniQPSolver<Matrix, Vector>::registry() {
        return QPSolverRegistry::instance();
    }

    template <class Matrix, class Vector>
    OmniQPSolver<Matrix, Vector>::OmniQPSolver() {
#ifdef UTOPIA_WITH_PETSC
        auto tron = utopia::make_unique<TaoQPSolver<Matrix, Vector>>();
        tron->tao_type(TAOTRON);
        tron->set_linear_solver(std::make_shared<GMRES<Matrix, Vector>>("bjacobi"));
        impl_ = std::move(tron);
#else
        impl_ = utopia::make_unique<SemismoothNewton<Matrix, Vector>>(
            QPSolverRegistry<Matrix, Vector>::default_linear_solver());
#endif
    }

    template <class Matrix, class Vector>
    OmniQPSolver<Matrix, Vector>::~OmniQPSolver() = default;

    template <class Matrix, class Vector>
    OmniQPSolver<Matrix, Vector> *OmniQPSolver<Matrix, Vector>::clone() const {
        auto cloned = utopia::make_unique<OmniQPSolver>();

        if (this->impl_) {
            cloned->impl_ = std::unique_ptr<QPSolver>(this->impl_->clone());
        }

        return cloned.release();
    }

    template <class Matrix, class Vector>
    bool OmniQPSolver<Matrix, Vector>::apply(const Vector &rhs, Vector &sol) {
        assert(static_cast<bool>(impl_));
        if (!impl_) return false;

        impl_->update(this->get_operator());
        impl_->set_box_constraints(this->get_box_constraints());
        return impl_->apply(rhs, sol);
    }

}  // namespace utopia

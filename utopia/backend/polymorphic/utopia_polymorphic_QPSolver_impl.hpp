#include "utopia_polymorphic_QPSolver.hpp"
#include "utopia_make_unique.hpp"
#include "utopia_polymorphic_LinearSolver.hpp"
#include "utopia_InputParameters.hpp"
#include "utopia_LinearSolverFactory.hpp"

// includes
// HOMEMADE
#include "utopia_SemismoothNewton.hpp"
#include "utopia_ProjectedGradient.hpp"
#include "utopia_ProjectedConjugateGradient.hpp"
#include "utopia_ProjectedGaussSeidel.hpp"

//PETSC
#ifdef WITH_PETSC
#include "utopia_petsc_Factorizations.hpp"
#include "utopia_petsc_TaoQPSolver.hpp"
#include "utopia_petsc_RowView.hpp"
#include "utopia_petsc_SemismoothNewton.hpp"
#include "petsctao.h"
#endif //WITH_PETSC



#include <map>
#include <string>
#include <functional>
#include <memory>

namespace utopia {


    //FIXME use factorymethod class
    template<class Matrix, class Vector>
    class QPSolverRegistry {
    public:
        using QPSolverPtr     = std::unique_ptr<utopia::QPSolver<Matrix, Vector>>;
        using LinearSolverPtr = std::unique_ptr<utopia::LinearSolver<Matrix, Vector>>;

        std::map<
            std::string,
            std::map<
                std::string,
                QPSolverPtr
            >
        > solvers;

        inline void register_solver(const std::string &backend,
                                    const std::string &type,
                                    QPSolverPtr &&solver
                                    )
        {
            solvers[backend][type] = std::move(solver);
        }

        QPSolverPtr find(const std::string &backend,
                         const std::string &type) const
        {
            auto b_it = solvers.find(backend);
            if(b_it == solvers.end()) {
                std::cerr << "[Error] backend " << backend << " not found using default" << std::endl;
                return default_solver();
            }

            auto s_it = b_it->second.find(type);

            if(s_it == b_it->second.end()) {
                std::cerr << "[Error] solver " << type << " not found using default" << std::endl;

                // if(b_it->second.empty()) {
                return default_solver();
                // } else {
                // 	return (*b_it->second.begin())->clone();
                // }
            }

            return QPSolverPtr(s_it->second->clone());
        }

        inline static LinearSolverPtr default_linear_solver()
        {
            auto ls = utopia::make_unique<PolymorphicLinearSolver<Matrix, Vector>>();
            InputParameters in;
            in.set("type", Solver::direct());
            ls->read(in);
            return std::move(ls);
// #ifdef WITH_PETSC
//             return utopia::make_unique<Factorization<Matrix, Vector>>();
// #else
//             return utopia::make_unique<BiCGStab<Matrix, Vector>>();
// #endif
        }

        inline static QPSolverPtr default_solver()
        {
            return utopia::make_unique<SemismoothNewton<Matrix, Vector>>(default_linear_solver());
        }

        inline static std::unique_ptr<QPSolverRegistry> make()
        {
            return utopia::make_unique<QPSolverRegistry>();
        }


        //FIXME make feature detection for registrations that always compile
        QPSolverRegistry()
        {
            register_solver("any", "ssnewton", utopia::make_unique<SemismoothNewton<Matrix, Vector>>(default_linear_solver()));
            register_solver("any", "pg",       utopia::make_unique<ProjectedGradient<Matrix, Vector>>());
            register_solver("any", "pcg",      utopia::make_unique<ProjectedConjugateGradient<Matrix, Vector>>());
            register_solver("any", "pgs",      utopia::make_unique<ProjectedGaussSeidel<Matrix, Vector>>());

#ifdef WITH_PETSC
            register_solver("petsc", "ssnewton", utopia::make_unique<SemismoothNewton<Matrix, Vector, PETSC_EXPERIMENTAL>>(default_linear_solver()));
            register_solver("petsc", "taoqp",    utopia::make_unique<TaoQPSolver<Matrix, Vector>>(default_linear_solver()));

            auto tron = utopia::make_unique<TaoQPSolver<Matrix, Vector>>();
            tron->tao_type(TAOTRON);
            tron->set_linear_solver(std::make_shared<GMRES<Matrix, Vector>>("bjacobi"));
            register_solver("petsc", TAOTRON, std::move(tron));
#endif //WITH_PETSC

        }

    };

    template<class Matrix, class Vector>
    void PolymorphicQPSolver<Matrix, Vector>::read(Input &in)
    {
        Super::read(in);

        std::string backend = Traits<Vector>::backend_info().get_name();
        std::string type  = "ssnewton";

        in.get("backend", backend);
        in.get("type",    type);

        auto registry = QPSolverRegistry<Matrix, Vector>::make();
        impl_ = registry->find(backend, type);
        impl_->read(in);
    }

    template<class Matrix, class Vector>
    PolymorphicQPSolver<Matrix, Vector>::PolymorphicQPSolver()
    {
#ifdef WITH_PETSC
        auto tron = utopia::make_unique<TaoQPSolver<Matrix, Vector>>();
        tron->tao_type(TAOTRON);
        tron->set_linear_solver(std::make_shared<GMRES<Matrix, Vector>>("bjacobi"));
        impl_ = std::move(tron);
#else
        impl_ = utopia::make_unique<SemismoothNewton<Matrix, Vector>>(QPSolverRegistry<Matrix, Vector>::default_linear_solver());
#endif

    }

    template<class Matrix, class Vector>
    PolymorphicQPSolver<Matrix, Vector>::~PolymorphicQPSolver()
    { }

    template<class Matrix, class Vector>
    PolymorphicQPSolver<Matrix, Vector> *  PolymorphicQPSolver<Matrix, Vector>::clone() const
    {
        auto cloned = utopia::make_unique<PolymorphicQPSolver>();

        if(this->impl_) {
            cloned->impl_ = std::unique_ptr<QPSolver<Matrix, Vector>>(this->impl_->clone());
        }

        return cloned.release();
    }

    template<class Matrix, class Vector>
    bool PolymorphicQPSolver<Matrix, Vector>::apply(const Vector &rhs, Vector &sol)
    {
        assert(static_cast<bool>(impl_));
        if(!impl_) return false;

        impl_->update(this->get_operator());
        impl_->set_box_constraints(this->get_box_constraints());
        return impl_->apply(rhs, sol);
    }

}



#ifndef UTOPIA_PETSC_TRUST_REGION_FACTORY_HPP
#define UTOPIA_PETSC_TRUST_REGION_FACTORY_HPP

#include "utopia_TrustRegionFactory.hpp"
#include "utopia_petsc_KSPTR.hpp"
#include "utopia_petsc_TaoSolver.hpp"
#include "utopia_FactoryMethod.hpp"

namespace utopia {
    /**
     * @brief      Front-end to create tr strategy objects.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<typename Matrix, typename Vector>
    class TRStrategyFactory<Matrix, Vector, PETSC> {
    public:
        typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblemT;
        typedef std::shared_ptr<TRSubproblemT> StrategyPtr;
        typedef utopia::IFactoryMethod<TRSubproblemT> FactoryMethodT;

        template<class Alg>
        using TRFactoryMethod = FactoryMethod<TRSubproblemT, Alg>;

        template<class Alg>
        using TRUnaryFactoryMethod = UnaryFactoryMethod<TRSubproblemT, std::string, Alg>;

        std::map<std::string, std::shared_ptr<FactoryMethodT> > strategies_;

        inline static StrategyPtr new_trust_region_strategy(const std::string & tag)
        {
            auto it = instance().strategies_.find(tag);
            if(it == instance().strategies_.end())
            {
                utopia_warning("Strategy not available, solving with SteihaugToint.");
                return std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();
            } else {
                return it->second->make();
            }
        }

    private:
        inline static const TRStrategyFactory &instance()
        {
            static TRStrategyFactory instance_;
            return instance_;
        }

        TRStrategyFactory()
        {
            init();
        }

        void init()
        {
            strategies_[Solver::cauchypoint()] 	  = std::make_shared< TRFactoryMethod< utopia::CauchyPoint<Matrix, Vector>> >();
            strategies_[Solver::dogleg()] 		  = std::make_shared< TRFactoryMethod< utopia::Dogleg<Matrix, Vector>> >();
            strategies_[Solver::steihaug_toint()] = std::make_shared< TRFactoryMethod< utopia::SteihaugToint<Matrix, Vector>> >();

            strategies_[Solver::automatic()] = std::make_shared< TRFactoryMethod< utopia::SteihaugToint<Matrix, Vector>> >();
            strategies_[Solver::nash()] 	 = std::make_shared< TRFactoryMethod< utopia::Nash<Matrix, Vector>> >();
            strategies_[Solver::lanczos()]   = std::make_shared< TRFactoryMethod< utopia::Lanczos<Matrix, Vector>> >();
            strategies_[Solver::cgne()]      = std::make_shared< TRFactoryMethod< utopia::CGNE<Matrix, Vector>> >();
        }
    };
}

#endif //UTOPIA_PETSC_TRUST_REGION_FACTORY_HPP

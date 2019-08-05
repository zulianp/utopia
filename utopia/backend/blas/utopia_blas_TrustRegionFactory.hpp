#ifndef UTOPIA_BLAS_TRUST_REGION_FACTORY_HPP
#define UTOPIA_BLAS_TRUST_REGION_FACTORY_HPP

#include "utopia_TrustRegionFactory.hpp"
#include "utopia_FactoryMethod.hpp"

namespace utopia {

    /**
     * @brief      Front-end to create tr strategy objects.
     *
     * @tparam     Matrix
     * @tparam     Vector
     */
    template<typename Matrix, typename Vector>
    class TRStrategyFactory<Matrix, Vector, BLAS> {
    public:
        typedef utopia::TRSubproblem<Matrix, Vector> TRSubproblemT;
        typedef std::shared_ptr<TRSubproblemT> StrategyPtr;
        typedef utopia::IFactoryMethod<TRSubproblemT> FactoryMethodT;

        template<class Alg>
        using TRFactoryMethod = FactoryMethod<TRSubproblemT, Alg>;
        std::map<std::string, std::shared_ptr<FactoryMethodT> > strategies_;

        inline static StrategyPtr new_trust_region_strategy(const std::string &tag)
        {
            auto it = instance().strategies_.find(tag);
            if(it == instance().strategies_.end()) {
                utopia_warning("Strategy not available, solving with SteihaugToint.");
                return std::make_shared< utopia::SteihaugToint<Matrix, Vector> >();
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
        }
    };
}

#endif //UTOPIA_BLAS_TRUST_REGION_FACTORY_HPP

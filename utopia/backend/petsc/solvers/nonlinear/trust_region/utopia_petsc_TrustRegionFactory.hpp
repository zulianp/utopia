#ifndef UTOPIA_PETSC_TRUST_REGION_FACTORY_HPP
#define UTOPIA_PETSC_TRUST_REGION_FACTORY_HPP 

#include "utopia_TrustRegionFactory.hpp"
#include "utopia_petsc_KSPTR.hpp"

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
		typedef std::shared_ptr< TRSubproblem<Matrix, Vector> > StrategyPtr;
		std::map<std::string, StrategyPtr> strategies_;

		inline static StrategyPtr new_trust_region_strategy(const std::string & tag)
		{
			auto it = instance().strategies_.find(tag);
			if(it == instance().strategies_.end()) 
			{
				std::cout<<"Strategy not available, solving with SteihaugToint.  \n";  
				return std::make_shared<utopia::SteihaugToint<Matrix, Vector> >(); 
			} 
			else 
			{
				return it->second;
			}
		}

	private:
		inline static const TRStrategyFactory &instance()
		{
			static TRStrategyFactory instance_;
			instance_.init();
			return instance_;
		}

		void init()
		{	
			strategies_[CAUCHYPOINT_TAG] 	= std::make_shared<utopia::CauchyPoint<Matrix, Vector> >(); 
			strategies_[DOGLEG_TAG] 		= std::make_shared<utopia::Dogleg<Matrix, Vector> >(); 
			strategies_[STEIHAUG_TOINT_TAG] = std::make_shared<utopia::KSP_TR<Matrix, Vector> >("stcg"); 
			strategies_[TOINT_TAG] 			= std::make_shared<utopia::KSP_TR<Matrix, Vector> >("qcg"); 
			strategies_[NASH_TAG] 			= std::make_shared<utopia::KSP_TR<Matrix, Vector> >("nash"); 
			strategies_[LANCZOS_TAG] 		= std::make_shared<utopia::KSP_TR<Matrix, Vector> >("gltr"); 
			strategies_[CGNE_TAG] 			= std::make_shared<utopia::KSP_TR<Matrix, Vector> >("cgne"); 
			strategies_[AUTO_TR_TAG] 		= std::make_shared<utopia::KSP_TR<Matrix, Vector> >(); 
		}
	};
}

#endif //UTOPIA_PETSC_TRUST_REGION_FACTORY_HPP

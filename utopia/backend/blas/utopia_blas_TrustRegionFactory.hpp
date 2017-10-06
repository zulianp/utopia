#ifndef UTOPIA_BLAS_TRUST_REGION_FACTORY_HPP
#define UTOPIA_BLAS_TRUST_REGION_FACTORY_HPP 

#include "utopia_TrustRegionFactory.hpp"

namespace utopia {

	/**
	 * @brief      Front-end to create tr strategy objects.
	 *
	 * @tparam     Matrix  
	 * @tparam     Vector  
	 */
	template<typename Matrix, typename Vector>
	class TRStrategyFactory<Matrix, Vector, BLAS>
	{
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
				strategies_[CAUCHYPOINT_TAG] 		= std::make_shared<utopia::CauchyPoint<Matrix, Vector> >(); 
				strategies_[DOGLEG_TAG] 			= std::make_shared<utopia::SteihaugToint<Matrix, Vector> >();  // THIS IS WRONG !!!!!!!
				strategies_[STEIHAUG_TOINT_TAG] 	= std::make_shared<utopia::SteihaugToint<Matrix, Vector> >(); 
			}
	};
}

#endif //UTOPIA_BLAS_TRUST_REGION_FACTORY_HPP

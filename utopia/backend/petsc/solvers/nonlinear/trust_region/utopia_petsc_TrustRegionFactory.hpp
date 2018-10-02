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
			strategies_[CAUCHYPOINT_TAG] 	= std::make_shared< TRFactoryMethod< utopia::CauchyPoint<Matrix, Vector>> >(); 
			strategies_[DOGLEG_TAG] 		= std::make_shared< TRFactoryMethod< utopia::Dogleg<Matrix, Vector>> >(); 
			strategies_[STEIHAUG_TOINT_TAG] = std::make_shared< TRFactoryMethod< utopia::SteihaugToint<Matrix, Vector>> >(); 


			strategies_[AUTO_TR_TAG] = std::make_shared< TRFactoryMethod< utopia::SteihaugToint<Matrix, Vector>> >(); 
			
			// strategies_[NASH_TAG] = std::make_shared< TRFactoryMethod< utopia::Nash<Matrix, Vector>> >(); 
			// strategies_[LANCZOS_TAG] = std::make_shared< TRFactoryMethod< utopia::Lanczos<Matrix, Vector>> >(); 
			// strategies_[CGNE_TAG] = std::make_shared< TRFactoryMethod< utopia::CGNE<Matrix, Vector>> >(); 


			// strategies_[AUTO_TR_TAG] 		= std::make_shared< TRFactoryMethod< utopia::SteihaugToint<Matrix, Vector>> >(); 
			
			// strategies_[TOINT_TAG] 			= std::make_shared< TRUnaryFactoryMethod< utopia::KSP_TR<Matrix, Vector>> >("qcg"); 
			// strategies_[NASH_TAG] 			= std::make_shared< TRUnaryFactoryMethod< utopia::KSP_TR<Matrix, Vector>> >("nash"); 
			// strategies_[LANCZOS_TAG] 		= std::make_shared< TRUnaryFactoryMethod< utopia::KSP_TR<Matrix, Vector>> >("gltr"); 
			// strategies_[CGNE_TAG] 			= std::make_shared< TRUnaryFactoryMethod< utopia::KSP_TR<Matrix, Vector>> >("cgne"); 
			
		}
	};
}

#endif //UTOPIA_PETSC_TRUST_REGION_FACTORY_HPP

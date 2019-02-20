#ifndef UTOPIA_EIKONAL_EQUATION_TEST_HPP
#define UTOPIA_EIKONAL_EQUATION_TEST_HPP 

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {
	class EikonalApp final : public FEApp {
	public:
		void run(Input &in) override;
		
		inline static std::string command()
		{
			return "-eikonal";
		}
	};
}


#endif //UTOPIA_EIKONAL_EQUATION_TEST_HPP

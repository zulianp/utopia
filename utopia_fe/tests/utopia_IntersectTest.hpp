#ifndef UTOPIA_INTERSECT_TEST_HPP
#define UTOPIA_INTERSECT_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

	class IntersectTest final : public FETest {
	public:
		void run(Input &in) override;

		inline static std::string command()
		{
			return "isect";
		}
	};

}

#endif //UTOPIA_INTERSECT_TEST_HPP

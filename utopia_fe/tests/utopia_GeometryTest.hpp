#ifndef UTOPIA_GEOMETRY_TEST_HPP
#define UTOPIA_GEOMETRY_TEST_HPP 

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

	class GeometryTest final : public FETest {
	public:
		void run(Input &in) override;

		inline static std::string command()
		{
			return "geo";
		}
	};

}

#endif //UTOPIA_GEOMETRY_TEST_HPP

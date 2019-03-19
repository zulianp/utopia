#ifndef UTOPIA_MECH_TEST_HPP
#define UTOPIA_MECH_TEST_HPP 

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

	class MechTest final : public FETest {
	public:
		void run(Input &in) override;

		inline static std::string command()
		{
			return "mech";
		}
	};

}

#endif //UTOPIA_MECH_TEST_HPP

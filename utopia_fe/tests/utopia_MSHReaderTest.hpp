#ifndef UTOPIA_MSH_READER_TEST_HPP
#define UTOPIA_MSH_READER_TEST_HPP 

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

	class MSHReaderTest final : public FETest {
	public:
		void run(Input &in) override;

		inline static std::string command()
		{
			return "msh";
		}
	};

}

#endif //UTOPIA_MSH_READER_TEST_HPP

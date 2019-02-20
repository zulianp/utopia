#ifndef UTOPIA_ASSEMBLY_TEST_HPP
#define UTOPIA_ASSEMBLY_TEST_HPP 

#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

	class AssemblyTest final : public FETest {
	public:
		void run(Input &in) override;

		inline static std::string command()
		{
			return "asm";
		}
	};

}

#endif //UTOPIA_ASSEMBLY_TEST_HPP

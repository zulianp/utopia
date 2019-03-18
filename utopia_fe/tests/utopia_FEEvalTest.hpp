#ifndef UTOPIA_FE_EVAL_TEST_HPP
#define UTOPIA_FE_EVAL_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

	class FEEvalTest final : public FETest {
	public:
		void run(Input &in) override;

		inline static std::string command()
		{
			return "fe_test";
		}
	};

}



#endif //UTOPIA_FE_EVAL_TEST_HPP

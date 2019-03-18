#ifndef UTOPIA_RMTR_APP_HPP
#define UTOPIA_RMTR_APP_HPP

#include <string>
#include "utopia_FEApp.hpp"

namespace utopia {
	class RMTRApp final : public FEApp {
	public:
		void run(Input &in) override;


		
		inline static std::string command()
		{
			return "-rmtr";
		}

		class SimulationInput;

	private:

		void solve_newton(const SimulationInput &in);
		void solve_rmtr(const SimulationInput &in);
	};
}

#endif //UTOPIA_RMTR_APP_HPP

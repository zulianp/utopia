#ifndef UTOPIA_APP_HPP
#define UTOPIA_APP_HPP

#include "utopia_ui.hpp"
#include <string>

namespace utopia {
	class App {
	public:
		virtual ~App() {}
		virtual void run(Input &in) = 0;
	};
}

#endif //UTOPIA_APP_HPP

#ifndef UTOPIA_APP_HPP
#define UTOPIA_APP_HPP

#include <string>

namespace utopia {
	class App {
	public:
		virtual ~App() {}
		virtual void run(const std::string &path) = 0;
	};
}

#endif //UTOPIA_APP_HPP

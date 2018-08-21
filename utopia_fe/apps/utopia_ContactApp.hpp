#ifndef UTOPIA_CONTACT_APP
#define UTOPIA_CONTACT_APP

#include "utopia_App.hpp"
#include <string>

namespace utopia {
	class ContactApp final : public App {
	public:
		void run(const std::string &path) override;
	};
}


#endif //UTOPIA_CONTACT_APP

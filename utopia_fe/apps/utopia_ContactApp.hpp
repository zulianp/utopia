#ifndef UTOPIA_CONTACT_APP
#define UTOPIA_CONTACT_APP

#include "utopia_LibMeshApp.hpp"
#include <string>

namespace utopia {

	class ContactApp final : public LibMeshApp {
	public:
		void run(const std::string &path) override;
	};
}


#endif //UTOPIA_CONTACT_APP

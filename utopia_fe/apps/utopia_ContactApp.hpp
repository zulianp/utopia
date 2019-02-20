#ifndef UTOPIA_CONTACT_APP
#define UTOPIA_CONTACT_APP

#include "utopia_FEApp.hpp"
#include "libmesh/parallel_mesh.h"
#include <string>


namespace utopia {

	class ContactApp final : public FEApp {
	public:
		void run(Input &in) override;

		inline static std::string command()
		{
			return "-contact";
		}
	};
}


#endif //UTOPIA_CONTACT_APP

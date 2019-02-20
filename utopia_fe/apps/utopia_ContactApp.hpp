#ifndef UTOPIA_CONTACT_APP
#define UTOPIA_CONTACT_APP

#include "utopia_FEApp.hpp"
#include "libmesh/parallel_mesh.h"
#include <string>


namespace utopia {

	class ContactApp final : public FEApp {
	public:
		void run(Input &in) override;
		void init(libMesh::Parallel::Communicator &comm) override;

		inline static std::string command()
		{
			return "-contact";
		}

	private:
		std::shared_ptr<libMesh::Parallel::Communicator> comm_;
	};
}


#endif //UTOPIA_CONTACT_APP

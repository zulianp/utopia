#ifndef UTOPIA_CONTACT_APP
#define UTOPIA_CONTACT_APP

#include "utopia_LibMeshApp.hpp"
#include "libmesh/parallel_mesh.h"
#include <string>


namespace utopia {

	class ContactApp final : public App {
	public:
		void run(const std::string &path) override;
		void init(libMesh::LibMeshInit &init);

		inline static std::string command()
		{
			return "-contact";
		}

	private:
		std::shared_ptr<libMesh::Parallel::Communicator> comm_;
	};
}


#endif //UTOPIA_CONTACT_APP

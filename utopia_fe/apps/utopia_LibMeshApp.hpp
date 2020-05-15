// #ifndef UTOPIA_LIBMESH_APP
// #define UTOPIA_LIBMESH_APP

// #include "utopia_App.hpp"
// #include "utopia_libmesh.hpp"

// #include <string>
// #include <memory>

// namespace utopia {

// 	class LibMeshApp : public App {
// 	public:
// 		virtual ~LibMeshApp() {}

// 		inline void init(const std::string &mesh_path)
// 		{
// 			comm_ = std::make_shared<libMesh::Parallel::Communicator>();
// 			mesh_ = std::make_shared<libMesh::DistributedMesh>(*comm_);
// 			mesh_->read(mesh_path);

// 			equation_systems_ = std::make_shared<libMesh::EquationSystems>(*mesh_);
// 		}

// 		inline std::shared_ptr<libMesh::Parallel::Communicator> comm()
// 		{
// 			return comm_;
// 		}

// 		inline std::shared_ptr<libMesh::DistributedMesh> mesh()
// 		{
// 			return mesh_;
// 		}

// 		inline std::shared_ptr<libMesh::EquationSystems> equation_systems()
// 		{
// 			return equation_systems_;
// 		}

// 	private:
// 		std::shared_ptr<libMesh::Parallel::Communicator> comm_;
// 		std::shared_ptr<libMesh::DistributedMesh> mesh_;
// 		std::shared_ptr<libMesh::EquationSystems> equation_systems_;
// 	};
// }

// #endif //UTOPIA_LIBMESH_APP

#ifndef UTOPIA_VOLUME_L2_PROJECTION_HPP
#define UTOPIA_VOLUME_L2_PROJECTION_HPP 

#include <memory>
#include <vector>
#include <utility>

#include "utopia.hpp"
#include "par_moonolith.hpp"

namespace libMesh {
	class MeshBase;
	class DofMap;
}

namespace moonolith {
	class Predicate;
	class MasterAndSlave;
}

namespace utopia {
	class FESpacesAdapter;

	// class FETransferAssembler {
	// public:

	// };

	class AssemblyData {
	public:
		virtual ~AssemblyData() {}
		//read
		//write
	};

	class FETransferParams {

	};

	class FETransfer {
	public:

		void init(
			const std::shared_ptr<libMesh::MeshBase> &master_mesh,
			const std::shared_ptr<libMesh::DofMap>   &master_dof_map,
			const std::shared_ptr<libMesh::MeshBase> &slave_mesh,
			const std::shared_ptr<libMesh::DofMap>   &slave_dof_map,
			const std::vector<std::pair<int, int>> &tags,
			const std::vector<std::pair<int, int>> &var_pairings = {{}}
			);

		FETransfer();
		~FETransfer();

		bool assemble();
		moonolith::SearchSettings search_settings;

	private:
		DSMatrixd coupling_;
		DSMatrixd mass_matrix_;
		DSMatrixd transfer_operator_;

		std::shared_ptr<libMesh::MeshBase> master_mesh;
		std::shared_ptr<libMesh::MeshBase> slave_mesh;
		
		std::shared_ptr<libMesh::DofMap> master_dof_map;
		std::shared_ptr<libMesh::DofMap> slave_dof_map;
		
		std::shared_ptr<moonolith::MasterAndSlave> predicate;
		std::shared_ptr<FESpacesAdapter> spaces;

		std::vector<std::pair<int, int>> var_pairings;

		//hacky stuff
		bool auto_tag;

		class Buffers;
		std::shared_ptr<Buffers> buffers_;

		template<int Dimension>
		bool assemble_aux();

		template<class TreeT>
		bool fill_tree(TreeT &tree);
	};
}

#endif //UTOPIA_VOLUME_L2_PROJECTION_HPP

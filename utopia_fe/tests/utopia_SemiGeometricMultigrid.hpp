#ifndef UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP
#define UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP 

#include "utopia_LinearSolver.hpp"
#include "utopia.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {
	class SemiGeometricMultigrid : public LinearSolver<DSMatrixd, DVectord> {
    public:
		void init(const LibMeshFunctionSpace &space, const std::size_t n_levels)
		{
			init(space.equation_systems(), n_levels);
		}

		void init(const libMesh::EquationSystems &es, const std::size_t n_levels);

		void update(const std::shared_ptr<const DSMatrixd> &op);
		bool apply(const DVectord &rhs, DVectord &sol);
		
		virtual void set_parameters(const Parameters params)
		{
			mg.set_parameters(params);
		}

		inline void verbose(const bool val)
		{
			mg.verbose(val);
		}

		SemiGeometricMultigrid();

	private:
		Multigrid<DSMatrixd, DVectord> mg;
		std::vector<std::shared_ptr<libMesh::UnstructuredMesh>> meshes;
		std::vector<std::shared_ptr<libMesh::EquationSystems>> equation_systems;
	};
}


#endif //UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP

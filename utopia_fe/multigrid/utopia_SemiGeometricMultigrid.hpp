#ifndef UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP
#define UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP 

#include "utopia_LinearSolver.hpp"
#include "utopia.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {
	class Contact;
	
	class SemiGeometricMultigrid : public LinearSolver<DSMatrixd, DVectord> {
    public:
    	typedef utopia::Multigrid<DSMatrixd, DVectord> MultigridT;
    	// typedef utopia::Multigrid<DSMatrixd, DVectord, PETSC_EXPERIMENTAL> MultigridT;

		void init(const LibMeshFunctionSpace &space, const std::size_t n_levels)
		{
			init(space.equation_systems(), n_levels);
		}

		void init(const libMesh::EquationSystems &es, const std::size_t n_levels);

		void update_contact(Contact &contact);
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

		inline void max_it(const unsigned int it) {
			mg.max_it(it);
		}

		void convert_to_block_solver()
		{
			is_block_solver_ = true;
		}

		SemiGeometricMultigrid(
			const std::shared_ptr<Smoother<DSMatrixd, DVectord> > &smoother = std::make_shared<GaussSeidel<DSMatrixd, DVectord>>(),
			const std::shared_ptr<LinearSolver<DSMatrixd, DVectord> > &linear_solver = std::make_shared<Factorization<DSMatrixd, DVectord>>()
		);

		inline MultigridT &algebraic()
		{
			return mg;
		}

	private:
		MultigridT mg;
		
		std::vector<std::shared_ptr<libMesh::UnstructuredMesh>> meshes;
		std::vector<std::shared_ptr<libMesh::EquationSystems>> equation_systems;

		std::vector<std::shared_ptr<DSMatrixd>> interpolators_;
		bool is_block_solver_;
	};
}


#endif //UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP

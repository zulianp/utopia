#ifndef UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP
#define UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP 

#include "utopia_libmesh_Types.hpp"
#include "utopia_IterativeSolver.hpp"
#include "utopia.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"

namespace utopia {
	class Contact;
	
	class SemiGeometricMultigrid final : public IterativeSolver<USparseMatrix, UVector> {
    public:
    	typedef utopia::Multigrid<USparseMatrix, UVector> MultigridT;
    	// typedef utopia::Multigrid<USparseMatrix, UVector, PETSC_EXPERIMENTAL> MultigridT;

		void init(const LibMeshFunctionSpace &space, const std::size_t n_levels)
		{
			init(space.equation_system(), n_levels);
		}

		void init(const libMesh::System &es, const std::size_t n_levels);

		void update_contact(Contact &contact);
		void update(const std::shared_ptr<const USparseMatrix> &op) override;
		bool apply(const UVector &rhs, UVector &sol) override;
		
		SemiGeometricMultigrid * clone() const override
		{
			return new SemiGeometricMultigrid();
		}

		virtual void set_parameters(const Parameters params) override
		{
			mg.set_parameters(params);
		}

		inline void verbose(const bool &val) override
		{
			mg.verbose(val);
		}

		inline void max_it(const SizeType &it) override {
			mg.max_it(it);
		}

		inline SizeType max_it() const override {
			return mg.max_it();
		}

		inline void atol(const double tol) {
			algebraic().atol(tol);
		}

		void convert_to_block_solver()
		{
			is_block_solver_ = true;
		}

		SemiGeometricMultigrid(
			const std::shared_ptr<Smoother<USparseMatrix, UVector> > &smoother = std::make_shared<GaussSeidel<USparseMatrix, UVector>>(),
			const std::shared_ptr<LinearSolver<USparseMatrix, UVector> > &linear_solver = std::make_shared<Factorization<USparseMatrix, UVector>>()
		);

		inline MultigridT &algebraic()
		{
			return mg;
		}

		void set_separate_subdomains(const bool val)
		{
			separate_subdomains_ = val;
		}

		void set_use_interpolation(const bool val)
		{
			use_interpolation_ = val;
		}

		void describe(std::ostream &os) const
		{
			os << "SemiGeometricMultigrid:\n";
			mg.describe(os);
		}

	private:
		MultigridT mg;
		
		std::vector<std::unique_ptr<libMesh::MeshBase>> meshes;
		std::vector<std::shared_ptr<libMesh::EquationSystems>> equation_systems;

		std::vector<std::shared_ptr<USparseMatrix>> interpolators_;
		std::vector<bool> use_interpolation_at_level_;
		bool is_block_solver_;
		bool separate_subdomains_;
		bool use_interpolation_;
		bool use_coarse_interpolators_;

		void generate_coarse_meshes(const libMesh::MeshBase &fine_mesh, const std::size_t n_levels, const int order_fine_level);
		std::unique_ptr<libMesh::MeshBase> generate_box_mesh(const libMesh::MeshBase &fine_mesh, const std::size_t n_levels);
	};
}


#endif //UTOPIA_FE_SEMI_GEOMETRIC_MULTIGRID_HPP

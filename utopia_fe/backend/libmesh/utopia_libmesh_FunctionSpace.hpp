#ifndef UTOPIA_LIBMESH_TAG_FUNCTION_SPACE_HPP
#define UTOPIA_LIBMESH_TAG_FUNCTION_SPACE_HPP

namespace utopia {

	class LibMeshFunctionSpace : public FunctionSpace<LibMeshFunctionSpace> {
	public:
		typedef libMesh::Real Scalar;
			
		inline LibMeshFunctionSpace(
			const std::shared_ptr<libMesh::EquationSystems> &equation_systems,
			const int system_num,
			const int subspace_id)
		: equation_systems_(equation_systems),
		  system_num_(system_num)
		{
			this->set_subspace_id(var_num);
		}
		
		inline libMesh::Order order()
		{
			return dof_map().variable_order(this->subspace_id());
		}

		inline libMesh::Order type()
		{
			return dof_map().variable_type(this->subspace_id());
		}
				
		inline libMesh::DofMap &dof_map() {
			return equation_system().get_dof_map();
		}
		
		inline const libMesh::DofMap &dof_map() const {
			return equation_system().get_dof_map();
		}

		inline libMesh::System &equation_system()
		{
			return equation_systems_->get_system(system_num_);
		}

		inline const libMesh::System &equation_system() const
		{
			return equation_systems_->get_system(system_num_);
		}
		
		inline libMesh::MeshBase &mesh() { return equation_systems_->get_mesh(); }
		inline const libMesh::MeshBase &mesh() const { return equation_systems_->get_mesh(); }
	
	private:
		std::shared_ptr<libMesh::EquationSystems> equation_systems_;
		int system_num_;
	};

	template<>
	class Traits<LibMeshFunctionSpace> {
	public:
		static const int Backend = LIBMESH_TAG;
		static const int Order = 1;
		static const int FILL_TYPE = FillType::DENSE;

		typedef double Scalar;
		typedef utopia::LMDenseMatrix Vector;
		typedef utopia::LMDenseVector Matrix;
		typedef libMesh::TensorValue<Scalar> TensorValueT;

		typedef utopia::LibMeshFunctionSpace Implementation;

		typedef libMesh::FEBase FE;

		typedef std::vector<std::vector<libMesh::FEBase::OutputShape>> FunctionType;
		typedef std::vector<std::vector<libMesh::FEBase::OutputGradient>> GradientType;
		typedef std::vector<std::vector<libMesh::FEBase::OutputDivergence>> DivergenceType;
		typedef std::vector<std::vector<TensorValueT>> JacobianType;
		typedef std::vector<std::vector<Vector>> CurlType;

		typedef std::vector<libMesh::Real> DXType;
	};

	typedef utopia::Traits<LibMeshFunctionSpace> LibMeshTraits;

	template<>
	class FunctionalTraits<LibMeshFunctionSpace, AssemblyContext<LIBMESH_TAG> > {
	public:
		inline static int type(const LibMeshFunctionSpace &space,  const AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return utopia::POLYNOMIAL_FUNCTION;
		}

		inline static int order(const LibMeshFunctionSpace &space, const AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return space.mesh().element_order(ctx.current_element);
		}
	};

}

#endif //UTOPIA_LIBMESH_TAG_FUNCTION_SPACE_HPP

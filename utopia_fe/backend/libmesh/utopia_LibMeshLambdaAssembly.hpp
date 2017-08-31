#ifndef UTOPIA_LIBMESH_LAMBDA_ASSEMBLY_HPP
#define UTOPIA_LIBMESH_LAMBDA_ASSEMBLY_HPP

#include "libmesh/linear_implicit_system.h"

namespace utopia {
	template<class Lambda>
	class LibMeshLambdaAssembly : public libMesh::LinearImplicitSystem::Assembly {
	public:
		LibMeshLambdaAssembly(const Lambda &lambda)
		: lambda_(lambda)
		{

		}
		virtual void assemble () 
		{
			lambda_();
		}

		Lambda lambda_;
	};

	template<class Lambda>
	LibMeshLambdaAssembly<Lambda> make_assembly(Lambda lambda)
	{
		return lambda;
	}
}

#endif //UTOPIA_LIBMESH_LAMBDA_ASSEMBLY_HPP

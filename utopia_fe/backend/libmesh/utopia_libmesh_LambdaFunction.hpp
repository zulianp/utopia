#ifndef UTOPIA_LIBMESH_LAMBDA_FUNCTION_HPP
#define UTOPIA_LIBMESH_LAMBDA_FUNCTION_HPP 

#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/analytic_function.h"

namespace utopia {

	template<typename Output>
	class LibMeshLambdaFunction : public libMesh::FunctionBase<Output> {
	public:
		typedef libMesh::Real Scalar;
		typedef libMesh::Point Point;
		
		typedef std::function<Scalar(const Point &)> ScalarFunction;
		typedef std::function<void(const Point &, libMesh::DenseVector<Output> &output)> VectorFunction;
		
		LibMeshLambdaFunction(const VectorFunction fun)
		: vector_fun_(fun)
		{}
		
		LibMeshLambdaFunction(const ScalarFunction fun)
		: scalar_fun_(fun)
		{}
		
		
		void init() override {}
		void clear() override {}
		
		libMesh::UniquePtr<libMesh::FunctionBase<Output> > clone() const override
		{
			//#ifdef LIBMESH_HAVE_CXX14_MAKE_UNIQUE
			//			using libMesh::make_unique;
			//			return make_unique< LibMeshLambdaFunction<Output> >(*this);
			//#else
			return libMesh::UniquePtr<libMesh::FunctionBase<Output> >( new LibMeshLambdaFunction<Output>(*this));
			//#endif
		}
		
		Output operator() (const Point &p, const Scalar time = 0.) override
		{
			return scalar_fun_(p);
		}
		
		void operator()(const Point &p, const Scalar time, libMesh::DenseVector<Output> &output) override
		{
			vector_fun_(p, output);
		}
		
	private:
		
		VectorFunction vector_fun_;
		ScalarFunction scalar_fun_;
	};
}


#endif //UTOPIA_LIBMESH_LAMBDA_FUNCTION_HPP

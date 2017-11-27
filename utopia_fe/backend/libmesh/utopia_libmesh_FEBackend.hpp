#ifndef UTOPIA_LIBMESH_FE_BACKEND_HPP
#define UTOPIA_LIBMESH_FE_BACKEND_HPP 

#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_AssemblyContext.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_LambdaFunction.hpp"

#include "utopia_FEBackend.hpp"
#include "utopia_FEBasicOverloads.hpp"


#include "libmesh/const_function.h"
#include <type_traits>

namespace utopia {

	template<typename T>
	inline static T inner(const libMesh::VectorValue<T> &left, const libMesh::VectorValue<T> &right)
	{
		return left * right;
	}

	template<typename T>
	inline static T inner(const libMesh::TensorValue<T> &left, const libMesh::TensorValue<T> &right)
	{
		return left.contract(right);
	}

	template<typename T>
	inline static T inner(const libMesh::DenseVector<T> &left, const libMesh::DenseVector<T> &right)
	{
		return left.dot(right);
	}



	template<typename T>
	inline static T inner(const libMesh::DenseVector<T> &left, const libMesh::VectorValue<T> &right)
	{
		T result = 0;
		std::size_t n = left.get_values().size();
		for(std::size_t i = 0; i < n; ++i) {
			result += left(i) * right(i);
		}

		return result;
	}

	template<typename T>
	inline static T inner(const libMesh::VectorValue<T> &left, const libMesh::DenseVector<T> &right)
	{
		T result = 0;
		std::size_t n = right.get_values().size();
		for(std::size_t i = 0; i < n; ++i) {
			result += right(i) * left(i);
		}

		return result;
	}

	template<typename T>
	inline static T inner(const libMesh::DenseMatrix<T> &left, const libMesh::DenseMatrix<T> &right)
	{
		std::size_t n = left.get_values().size();

		T result = 0;
		for(std::size_t i = 0; i < n; ++i) {
			result += left.get_values()[i] * right.get_values()[i];
		}

		return result;
	}

	template<class Left, typename T>
	inline static T inner(const Wrapper<Left, 1> &left, const libMesh::VectorValue<T> &right)
	{
		return inner(left.implementation(), right);
	}

	template<class Left, typename T>
	inline static T inner(const Wrapper<Left, 1> &left, const libMesh::DenseVector<T> &right)
	{
		return inner(left.implementation(), right);
	}


	template<typename T, class Right>
	inline static T inner(const libMesh::VectorValue<T> &left, const  Wrapper<Right, 1> &right)
	{
		return inner(right.implementation(), left);
	}

	template<typename T, class Right>
	inline static T inner(const libMesh::DenseVector<T> &left, const Wrapper<Right, 1> &right)
	{
		return inner(right.implementation(), left);
	}


	template<>
	class FEBackend<LIBMESH_TAG> {
	public:
		typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;

		static const int Backend = LIBMESH_TAG;
		static const int Order = 1;
		static const int FILL_TYPE = FillType::DENSE;

		typedef TraitsT::Scalar Scalar;
		typedef TraitsT::Vector Vector;
		typedef TraitsT::Matrix Matrix;
		typedef TraitsT::TensorValueT TensorValueT;
		typedef TraitsT::VectorValueT VectorValueT;
		typedef TraitsT::FE FE;
		typedef TraitsT::FunctionType FunctionType;
		typedef TraitsT::GradientType GradientType;
		typedef TraitsT::DivergenceType DivergenceType;
		typedef TraitsT::JacobianType JacobianType;
		typedef TraitsT::CurlType CurlType;
		typedef TraitsT::DXType DXType;

		typedef libMesh::DenseVector<libMesh::Real> DenseVectorT;
		typedef std::vector< std::vector<DenseVectorT> > VectorFunctionType;

		//////////////////////////////////////////////////////////////////////////////////////////
		
		//system of equations
		template<class... Eq>
		static void init_context(const Equations<Eq...> &eqs, AssemblyContext<LIBMESH_TAG> &ctx)
		{

		}

		//single equations
		template<class Left, class Right>
		static void init_context(const Equality<Left, Right> &eq, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			
		}

		//general expressions that might contain multilinear forms
		template<class Expr>
		static void init_context(const Expression<Expr> &expr, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			
		}

		template<class Tensor, int Order>
		static void init_tensor(Wrapper<Tensor, Order> &t, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			ctx.init_tensor(t, true);
		}


		class ConstraintsInitializer {
		public:

			template<class Constr>
			inline void operator()(const int, const Constr &constr) const
			{
				std::cout << "unimplemented constraint exprassion" << std::endl;
				assert(false);
			}

			template<class F, typename T>
			inline void operator()(const int, const DirichletBoundaryCondition<
																Equality<
																	TrialFunction<LibMeshFunctionSpace>,
											 						FunctionCoefficient<F, T, 0>
											 						>
											 					> &cond) const
			{
				std::vector<unsigned int> vars(1); vars[0] = cond.expr().left().space_ptr()->subspace_id();
				std::set<libMesh::boundary_id_type> bt;
				bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
				libMesh::DirichletBoundary d_bc(bt, vars, LibMeshLambdaFunction<libMesh::Real>(cond.expr().right().fun()) );
				cond.expr().left().space_ptr()->dof_map().add_dirichlet_boundary(d_bc);
			}

			template<class T>
			inline void operator()(const int, const DirichletBoundaryCondition<
																Equality<TrialFunction<LibMeshFunctionSpace>,
																		ConstantCoefficient<T, 0> >
														> &&cond)
			{
				std::vector<unsigned int> vars(1); vars[0] = cond.expr().left().space_ptr()->subspace_id();
				std::set<libMesh::boundary_id_type> bt;
				bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
				libMesh::DirichletBoundary d_bc(bt, vars, libMesh::ConstFunction<libMesh::Real>(cond.expr().right()) );
				cond.expr().left().space_ptr()->dof_map().add_dirichlet_boundary(d_bc);
			}

			template<class T>
			void strong_enforce( DirichletBoundaryCondition<
								Equality<TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>>,
								ConstantCoefficient<T, 1> >
								> &&cond)
			{
				std::cout << "unimplemented constraint exprassion" << std::endl;
				assert(false);

				// std::vector<unsigned int> vars(1); vars[0] = cond.expr().left().space_ptr()->subspace_id();
				// std::set<libMesh::boundary_id_type> bt;
				// bt.insert(cond.boundary_tags().begin(), cond.boundary_tags().end());
				// libMesh::DirichletBoundary d_bc(bt, vars, libMesh::ConstFunction<T>(cond.expr().right()) );
				// cond.expr().left().space_ptr()->dof_map().add_dirichlet_boundary(d_bc);
			}

		};

		template<class... Constr>
		static void init_constraints(const FEConstraints<Constr...> &constr)
		{

			ConstraintsInitializer c_init;
			constr.each(c_init);
			/*
				template<class Right, typename T>
				void strong_enforce(const DirichletBoundaryCondition<
									Equality<LibMeshFEFunction,
									FunctionCoefficient<Right, T, 0> >
									> &cond)
				{
				
				}

			*/

		}


		template<typename T>
		inline static const T &get(const std::vector<std::vector<T> > &v, const std::size_t qp, const std::size_t i)
		{
			return v[i][qp];
		}

		template<typename T>
		inline static const T &get(const std::vector<T> &v, const std::size_t qp, const std::size_t)
		{
			return v[qp];
		}

		template<typename T, int Order>
		inline static const T get(const ConstantCoefficient<T, Order> &c, const std::size_t qp, const std::size_t i)
		{
			return c[i];
		}

		template<typename T>
		inline static void add(Matrix &mat, const int i, const int j, const T value)
		{
			mat.add(i, j, value);
		}

		template<typename T>
		inline static void add(Vector &vec, const int i, const int j, const T value)
		{
			assert(j == 0);
			vec.add(i, value);
		}

		//////////////////////////////////////////////////////////////////////////////////////////
		static inline unsigned int offset(const LibMeshFunctionSpace &space, const AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return ctx.offset[space.subspace_id()];
		}


		static inline unsigned int n_shape_functions(const LibMeshFunctionSpace &space, const AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return ctx.fe()[space.subspace_id()]->n_shape_functions();
		}

		static inline unsigned int n_shape_functions(const ProductFunctionSpace<LibMeshFunctionSpace> &space, const AssemblyContext<LIBMESH_TAG> &ctx)
		{
			unsigned int ret = 0;

			space.each([&](const int, const LibMeshFunctionSpace &subspace) {
				ret += ctx.fe()[subspace.subspace_id()]->n_shape_functions();
			});

			return ret;
		}
		//////////////////////////////////////////////////////////////////////////////////////////

		//Function
		// scalar fe functions
		static const FunctionType &fun(const TrialFunction<LibMeshFunctionSpace> &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			assert(!ctx.trial().empty());
			return ctx.trial()[fun.space_ptr()->subspace_id()]->get_phi();
		}

		static const FunctionType &fun(const TestFunction<LibMeshFunctionSpace> &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			assert(!ctx.test().empty());
			return ctx.test()[fun.space_ptr()->subspace_id()]->get_phi();
		}

		// vector fe functions
		static void fun_aux(
			ProductFunctionSpace<LibMeshFunctionSpace> &space,
			std::vector<std::unique_ptr<FE> > &fe_object,
			AssemblyContext<LIBMESH_TAG> &ctx,
			VectorFunctionType &ret)
		{
			assert(!fe_object.empty());

			const auto &sub_0 = space[0];

			const std::size_t n_quad_points = fe_object[sub_0.subspace_id()]->get_phi()[0].size();

			unsigned int n_shape_functions = 0;
			space.each([&n_shape_functions, &fe_object](const int index, const LibMeshFunctionSpace &subspace) {
				n_shape_functions += fe_object[subspace.subspace_id()]->n_shape_functions();
			});

			ret.resize(n_shape_functions); 
			for(std::size_t i = 0; i < n_shape_functions; ++i) {
				ret[i].resize(n_quad_points);
				for(auto &r : ret[i]) {
					r.resize(space.n_subspaces());
					r.zero();
				}
			}

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				int offset = 0;
				space.each([&offset, &fe_object, &ret, qp](const int sub_index, const LibMeshFunctionSpace &s) {
					const auto &fe = fe_object[s.subspace_id()];
					assert((static_cast<bool>(fe)));
					const auto & fun = fe->get_phi();
					uint n_shape_i = fe->n_shape_functions();

					for(uint j = 0; j < n_shape_i; ++j) {
						ret[offset++][qp](sub_index) = fun[j][qp];
					}
				});
			}
		}

		static VectorFunctionType fun(const TestFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			VectorFunctionType ret;

			fun_aux(*space_ptr, ctx.test(), ctx, ret);
			return ret;
		}

		static VectorFunctionType fun(const TrialFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			VectorFunctionType ret;
			fun_aux(*space_ptr, ctx.trial(), ctx, ret);
			return ret;
		}

		//////////////////////////////////////////////////////////////////////////////////////////

		//Gradient
		// scalar fe functions
		static const GradientType &grad(const TrialFunction<LibMeshFunctionSpace> &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return ctx.trial()[fun.space_ptr()->subspace_id()]->get_dphi();
		}

		static const GradientType &grad(const TestFunction<LibMeshFunctionSpace> &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return ctx.test()[fun.space_ptr()->subspace_id()]->get_dphi();
		}

		// vector fe functions
		static void grad_aux(
			ProductFunctionSpace<LibMeshFunctionSpace> &space,
			std::vector<std::unique_ptr<FE> > &fe_object,
			AssemblyContext<LIBMESH_TAG> &ctx,
			JacobianType &ret)
		{

			assert(!fe_object.empty());
			const auto &sub_0 = space[0];
			const std::size_t n_quad_points = fe_object[sub_0.subspace_id()]->get_phi()[0].size();
			const uint dim = space[0].mesh().mesh_dimension();

			uint n_shape_functions = 0;
			space.each([&fe_object, &n_shape_functions](const int, const LibMeshFunctionSpace &subspace) {
				n_shape_functions += fe_object[subspace.subspace_id()]->n_shape_functions();
			});

			ret.resize(n_shape_functions); 
			for(std::size_t i = 0; i < n_shape_functions; ++i) {
				ret[i].resize(n_quad_points);
				//TensorValue is by default initialized to 0s
			}

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				int offset = 0;
				space.each([&offset, &fe_object, &ret, &space, qp, dim](const int sub_index, const LibMeshFunctionSpace &s) {
					const auto &fe = fe_object[s.subspace_id()];	
					const uint n_shape_i = fe->n_shape_functions();
					
					for(uint j = 0; j < n_shape_i; ++j, offset++) {
						const auto &grad = fe->get_dphi()[j][qp];

						for(uint d = 0; d < dim; ++d) {
							ret[offset][qp](sub_index, d) = grad(d);
						}
					}
				});
			}
		}

		static JacobianType grad(const TrialFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			assert(!ctx.trial().empty());

			auto space_ptr = fun.space_ptr();
			JacobianType ret;
			grad_aux(*space_ptr, ctx.trial(), ctx, ret);
			return ret;
		}

		static JacobianType grad(const TestFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			assert(!ctx.test().empty());

			auto space_ptr = fun.space_ptr();
			JacobianType ret;
			grad_aux(*space_ptr, ctx.test(), ctx, ret);
			return ret;
		}


			//////////////////////////////////////////////////////////////////////////////////////////

			//Divergence
			// vector fe functions
		static void div_aux(
			ProductFunctionSpace<LibMeshFunctionSpace> &space,
			std::vector<std::unique_ptr<FE> > &fe_object,
			AssemblyContext<LIBMESH_TAG> &ctx,
			FunctionType &ret)
		{
			const auto &sub_0 = space[0];
			ret.resize(fe_object[sub_0.subspace_id()]->get_phi().size()); 

			assert(!fe_object.empty());
			const std::size_t n_quad_points = fe_object[sub_0.subspace_id()]->get_phi()[0].size();
			const uint dim = space[0].mesh().mesh_dimension();


			uint n_shape_functions = 0;
			space.each([&fe_object, &n_shape_functions](const int, const LibMeshFunctionSpace &subspace) {
				n_shape_functions += fe_object[subspace.subspace_id()]->n_shape_functions();
			});

			ret.resize(n_shape_functions); 
			for(std::size_t i = 0; i < n_shape_functions; ++i) {
				ret[i].resize(n_quad_points);
				//TensorValue is by default initialized to 0s
			}

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				int offset = 0;
				space.each([&offset, &fe_object, &ret, &space, qp, dim](const int sub_index, const LibMeshFunctionSpace &s) {
					const auto &fe = fe_object[s.subspace_id()];	
					const uint n_shape_i = fe->n_shape_functions();
					
					for(uint j = 0; j < n_shape_i; ++j, offset++) {
						const auto &grad = fe->get_dphi()[j][qp];
						ret[offset][qp] = grad(sub_index);
					}
				});
			}
		}

		static FunctionType div(const TrialFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			FunctionType ret;
			div_aux(*space_ptr, ctx.trial(), ctx, ret);
			return ret;
		}

		static FunctionType div(const TestFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			FunctionType ret;
			div_aux(*space_ptr, ctx.test(), ctx, ret);
			return ret;
		}


			//////////////////////////////////////////////////////////////////////////////////////////
			//Curl
		template<class CurlT>
		static void eval_curl(const uint dim, const TensorValueT &deriv, CurlT &result)
		{	
			switch(dim) {
				case 2:
				{
					result(2) = deriv(1, 0) - deriv(0, 1);
					break;
				}

				case 3:
				{
					result(0) = deriv(2, 1) - deriv(1, 2);
					result(1) = deriv(0, 2) - deriv(2, 0);
					result(2) = deriv(1, 0) - deriv(0, 1);
					break;
				}

				default :
				{
					assert(false);
				}
			}
		}

		static void curl_aux(
			ProductFunctionSpace<LibMeshFunctionSpace> &space,
			std::vector<std::unique_ptr<FE> > &fe_object,
			AssemblyContext<LIBMESH_TAG> &ctx,
			CurlType &ret)
		{
				const auto &sub_0 = space[0];
				JacobianType grads;

				grad_aux(
					space,
					fe_object,
					ctx,
					grads);

				const uint n_shape_x = fe_object[sub_0.subspace_id()]->n_shape_functions();
				uint n_shape_functions = 0;

				space.each([&n_shape_functions, &fe_object, &n_shape_x](const int, LibMeshFunctionSpace &subspace){
					assert(n_shape_x == fe_object[subspace.subspace_id()]->n_shape_functions());
					n_shape_functions += fe_object[subspace.subspace_id()]->n_shape_functions();
				});

				const std::size_t n_subspaces = space.n_subspaces();
				const std::size_t n_quad_points = grads[0].size();
				const std::size_t dim = sub_0.mesh().mesh_dimension();

				ret.resize(n_shape_functions); 
				for(auto &r : ret) {
					r.resize(n_quad_points);
				}

				for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
					uint offset = 0;
					for(uint i = 0; i < n_subspaces; ++i) {
						for(uint j = 0; j < n_shape_x; ++j, ++offset) {
							auto &ret_j  = ret[offset][qp];
							auto &grad_j = grads[offset][qp];
							eval_curl(dim, grad_j, ret_j);
						}
					}
				}
		}


		static CurlType curl(const TrialFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			CurlType ret;
			curl_aux(*space_ptr, ctx.trial(), ctx, ret);
			return ret;
		}

		static CurlType curl(const TestFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			CurlType ret;
			curl_aux(*space_ptr, ctx.test(), ctx, ret);
			return ret;
		}


		//////////////////////////////////////////////////////////////////////////////////////////
		//Interpolate
		template<class Tensor>
		static void gather_interp_values(
			const Interpolate<Wrapper<Tensor, 1>, TrialFunction<LibMeshFunctionSpace> > &interp,
			Vector &element_values,
			AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto &c   = interp.coefficient();
			auto &f   = interp.fun();

			auto space_ptr = f.space_ptr();
			const auto &mesh = space_ptr->mesh();
			const auto &elem_ptr = mesh.elem(ctx.current_element());
			const auto &dof_map = space_ptr->dof_map();

			std::vector<libMesh::dof_id_type> indices;
			dof_map.dof_indices(elem_ptr, indices, space_ptr->subspace_id());

			element_values = zeros(indices.size());

			Write<Vector> w(element_values);
			Read<Wrapper<Tensor, 1>> r(c);

			for(std::size_t i = 0; i < indices.size(); ++i) {
				element_values.set(i, c.get(indices[i]));
			}
		}

		template<class Tensor, class Space>
		static std::vector<Scalar> fun(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<Space> > &interp, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto &c  = interp.coefficient();
			auto &f  = interp.fun();
			auto &&g = fun(f, ctx);

			Vector element_values;
			gather_interp_values(interp, element_values, ctx);

			const std::size_t n_shape_functions = g.size();
			const std::size_t n_quad_points = g[0].size();

			std::vector<Scalar> ret(n_quad_points, 0.);

			for(std::size_t i = 0; i < n_shape_functions; ++i) {
				for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
					ret[qp] += g[i][qp] * element_values.get(i);
				}
			}

			return ret;
		}
	
		template<class Tensor, class Space>
		static std::vector<DenseVectorT> fun(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<ProductFunctionSpace<Space>> > &interp, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto &c  = interp.coefficient();
			auto &f  = interp.fun();
			auto &&g = fun(f, ctx);

			Vector element_values;
			gather_interp_values(interp, element_values, ctx);

			const std::size_t n_shape_functions = g.size();
			const std::size_t n_quad_points = g[0].size();

			std::vector<DenseVectorT> ret(n_quad_points);

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				ret[qp].resize(f.codim());
				ret[qp].zero();

				for(std::size_t i = 0; i < n_shape_functions; ++i) {
					if(std::is_rvalue_reference<decltype(g)>::value) {
						g[i][qp] *= element_values.get(i);
						ret[qp] += g[i][qp];
					} else {
						auto temp = g[i][qp];
						temp *= element_values.get(i);
						ret[qp] += temp;
					}
				}
			}

			return ret;
		}

		template<class Tensor>
		static std::vector<Vector> grad(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<LibMeshFunctionSpace> > &interp, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto &c  = interp.coefficient();
			auto &f  = interp.fun();
			auto &&g = grad(f, ctx);

			Vector element_values;
			gather_interp_values(interp, element_values, ctx);

			const std::size_t n_shape_functions = g.size();
			const std::size_t n_quad_points = g[0].size();

			std::vector<Vector> ret(n_quad_points);

			const uint dim = f.space_ptr()->mesh().mesh_dimension();
			
			for(auto &r : ret) {
				r = zeros(dim);
			}

			for(std::size_t i = 0; i < n_shape_functions; ++i) {
				for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
					// ret[qp] += element_values.get(i) * g[i][qp];
					add(ret[qp], element_values.get(i) * g[i][qp]);
				}
			}

			return ret;
		}

		static void add(Matrix &mat, const TensorValueT &t)
		{
			auto s = size(mat);
			Write<Matrix> w_m(mat);

			for(int i = 0; i < s.get(0); ++i) {
				for(int j = 0; j < s.get(1); ++j) {
					mat.add(i, j, t(i, j));
				}
			}
		}

		static void add(Vector &v, const VectorValueT &t)
		{
			auto s = size(v);
			Write<Vector> w_v(v);

			for(int i = 0; i < s.get(0); ++i) {
					v.add(i, t(i));
			}
		}


		template<class Tensor>
		static std::vector<Matrix> grad(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace> > > &interp, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto &c   = interp.coefficient();
			auto &f   = interp.fun();
			auto space_ptr = f.space_ptr();

			auto &&g  = grad(f, ctx);



			Vector element_values;
			gather_interp_values(interp, element_values, ctx);

			const SizeType rows = space_ptr->n_subspaces();
			const SizeType cols = space_ptr->subspace(0).mesh().mesh_dimension();
			Size s{rows, cols};

			
			const std::size_t n_shape_functions = g.size();
			const std::size_t n_quad_points = g[0].size();

			std::vector<Matrix> ret(n_quad_points);
			
			for(auto &r : ret) {
				r = zeros(s);
			}

			for(std::size_t i = 0; i < n_shape_functions; ++i) {
				for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
					add(ret[qp], element_values.get(i) * g[i][qp]);
				}
			}

			return ret;
		}

		template<class Tensor>
		static void gather_interp_values(
			const Interpolate<Wrapper<Tensor, 1>, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> > &interp,
			Vector &element_values,
			AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto &c   = interp.coefficient();
			auto &f   = interp.fun();

			auto space_ptr = f.space_ptr();
			const auto &sub_0 = space_ptr->subspace(0);
			const auto &mesh = sub_0.mesh();
			const auto &dof_map  = sub_0.dof_map();

			const auto &elem_ptr = mesh.elem(ctx.current_element());
			
			std::vector<libMesh::dof_id_type> prod_indices;
			std::vector<libMesh::dof_id_type> indices;

			space_ptr->each([&](const int sub_index, const LibMeshFunctionSpace &space) {
				dof_map.dof_indices(elem_ptr, indices, space_ptr->subspace_id());
				prod_indices.insert(prod_indices.end(), indices.begin(), indices.end());
			});

			const std::size_t n_indices = prod_indices.size();
			element_values = zeros(n_indices);

			Write<Vector> w(element_values);
			Read<Wrapper<Tensor, 1>> r(c);

			for(std::size_t i = 0; i < n_indices; ++i) {
				element_values.set(i, c.get(prod_indices[i]));
			}
		}

		//////////////////////////////////////////////////////////////////////////////////////////

		inline static auto apply_binary(const Number<double> &left, const Gradient<TrialFunction<LibMeshFunctionSpace> > &right, const Multiplies &op, AssemblyContext<LIBMESH_TAG> &ctx) -> GradientType
		{
			return multiply(static_cast<double>(left), right, ctx);
		}

		template<class Space, class Op>
		inline static auto apply_binary(const Number<double> &left, const TrialFunction<Space> &right, const Op &op, AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(fun(right, ctx))>::type 
		{
			typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			const double num = left;
			for(auto &v : ret) {
				for(auto &s : v) {
					s = op.apply(num, s);
				}
			}

			return ret;
		}

		template<class Space, class Op>
		inline static auto apply_binary(const Number<double> &left, const TestFunction<Space> &right, const Op &op, AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(fun(right, ctx))>::type 
		{
			typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			const double num = left;
			for(auto &v : ret) {
				for(auto &s : v) {
					s = op.apply(num, s);
				}
			}

			return ret;
		}

		//alpha * (grad_t + grad)
		template<class Left, class Right>
		inline static JacobianType grad_t_plus_grad(const double scaling, const Left &left, const Right &right, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			JacobianType ret  = grad(right, ctx);
			auto && grad_left = grad(left,  ctx);

			for(std::size_t i = 0; i < ret.size(); ++i) {
				for(std::size_t qp = 0; qp < ret[i].size(); ++qp) {
					ret[i][qp] += grad_left[i][qp].transpose();
					ret[i][qp] *= scaling;
				}
			}

			return ret;
		}

			//////////////////////////////////////////////////////////////////////////////////////////////////////

		inline static auto multiply(
			const std::vector<std::vector<double>> &left,
			const std::vector<Scalar> &right,
			const AssemblyContext<LIBMESH_TAG> &ctx
			) -> std::vector<std::vector<double>> 
		{
			auto ret = left;
			for(auto &r : ret) {
				for(std::size_t i = 0; i < right.size(); ++i) {
					r[i] *= right[i];
				}
			}

			return ret;
		}

		template<class Space>
		inline static auto multiply(
			const Scalar &left,
			const Gradient<TrialFunction<Space> > &right, 
			AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type 
		{
			typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			for(auto &v : ret) {
				for(auto &s : v) {
					s = left * s;
				}
			}

			return ret;
		}

		// template<class Space>
		inline static auto multiply(
			const LMDenseMatrix &left,
			const Gradient<TrialFunction<LibMeshFunctionSpace> > &right, 
			AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type 
		{
			typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			Read<LMDenseMatrix> r_l(left);

			Size s_l = size(left);

			for(auto &v : ret) {
				for(auto &s : v) {
					auto s_copy = s;

					for(uint i = 0; i < s_l.get(0); ++i) {
						s(i) = left.get(i, 0) * s_copy(0);

						for(uint j = 1; j < s_l.get(1); ++j) {
							s(i) += left.get(i, j) * s_copy(j);
						}	
					}

				}
			}

			return ret;
		}

		inline static void multiply(const LMDenseMatrix &left, const TensorValueT &right, TensorValueT &out)
		{
			Size s_l = size(left);
			Read<LMDenseMatrix> r_l(left);

			for(uint k = 0; k < LIBMESH_DIM; ++k) {
				for(uint i = 0; i < s_l.get(0); ++i) {
					out(i, k) = left.get(i, 0) * right(0, k);

					for(uint j = 1; j < s_l.get(1); ++j) {
						out(i, k) += left.get(i, j) * right(j, k);
					}	
				}
			}
		}

		inline static void multiply(const LMDenseMatrix &left, const DenseVectorT &right, DenseVectorT &out)
		{
			left.implementation().vector_mult(out, right);
		}

		// template<class Space>
		inline static auto multiply(
			const LMDenseMatrix &left,
			const Gradient<TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> > &right, 
			AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type 
		{
			typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			

			

			for(auto &v : ret) {
				for(auto &s : v) {
					auto s_copy = s;
					multiply(left, s_copy, s);
				}
			}

			return ret;
		}

		template<class Space>
		inline static auto multiply(
			const std::vector<double> &left,
			const Gradient<TrialFunction<Space> > &right, 
			AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type 
		{
			typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			const std::size_t n_functions = ret.size();
			const std::size_t n_quad_points = ret[0].size();

			assert(n_quad_points == left.size());

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					ret[i][qp] *= left[qp];
				}
			}

			return ret;
		}

		template<class Left, class Space>
		inline static auto multiply(
			const Left &left,
			const Gradient<TestFunction<Space> > &right, 
			AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type 
		{
			typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			const std::size_t n_functions = right.size();
			const std::size_t n_quad_points = right[0].size();

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					ret[i][qp] = left * ret[i][qp];
				}
			}


			return ret;
		}

		template<class Space>
		inline static auto multiply(const std::vector<double> &left, const TrialFunction<Space> &right,  AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(fun(right, ctx))>::type
		{
			typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			const std::size_t n_quad_points = left.size();
			const std::size_t n_functions = right.size();

			assert(n_quad_points == left.size());

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					ret[i][qp] *= left[qp];
				}
			}

			return ret;
		}

		template<class Space>
		inline static auto multiply(const std::vector<Matrix> &left, const Gradient<TrialFunction<Space> > &right,  AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type
		{
			typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			const std::size_t n_quad_points = left.size();
			const std::size_t n_functions = ret.size();

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					auto g = ret[i][qp];
					multiply(left[qp], g, ret[i][qp]);
				}
			}

			return ret;
		}

		template<typename T1, typename T2>
		inline static auto multiply(const std::vector<std::vector<T1>> &left, const std::vector<std::vector<T2>> &right, AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<std::vector<decltype(T1() * T2())>>
		{
			std::vector<std::vector<decltype(T1() * T2())>> ret(left.size());

			return ret;
		}

		template<class Space>
		inline static auto multiply(const std::vector<Matrix> &left, const TrialFunction<Space> &right, AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(fun(right, ctx))>::type
		{
			typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			const std::size_t n_quad_points = left.size();
			const std::size_t n_functions = ret.size();

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					auto v = ret[i][qp];
					multiply(left[qp], v, ret[i][qp]);
				}
			}

			return ret;
		}


		inline static auto multiply(const std::vector<Vector> &left, const TrialFunction<LibMeshFunctionSpace> &right, AssemblyContext<LIBMESH_TAG> &ctx) -> std::vector<std::vector<Vector>> 
		{
			const auto &f = fun(right, ctx);
			std::vector<std::vector<Vector>> ret(f.size());

			const std::size_t n_quad_points = left.size();
			const std::size_t n_functions = ret.size();
			
			for(std::size_t i = 0; i < n_functions; ++i) {
				ret[i].resize(left.size());
				for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
						// multiply(left[qp], f[i][qp], ret[i][qp]);
					ret[i][qp] = f[i][qp] * left[qp];
				}
			}

			return ret;
		}

		template<class Left, class Space>
		inline static auto multiply(const Left &left, const TrialFunction<Space> &right, AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(fun(right, ctx))>::type
		{
			typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			for(auto &v : ret) {
				for(auto &s : v) {
					s = left * s;
				}
			}

			return ret;
		}

		template<class Left, class Space>
		inline static auto multiply(const Number<Left> &left, const TestFunction<Space> &right,  AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(fun(right, ctx))>::type
		{
			typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			const std::size_t n_functions = ret.size();
			const std::size_t n_quad_points = ret[0].size();

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					ret[i][qp] = left * ret[i][qp];
				}
			}

			return ret;
		}

		template<class Space>
		inline static auto multiply(const Matrix &left, const TestFunction<ProductFunctionSpace<Space>> &right,  AssemblyContext<LIBMESH_TAG> &ctx) -> VectorFunctionType
		{
			VectorFunctionType ret = fun(right, ctx);

			const std::size_t n_functions = ret.size();
			const std::size_t n_quad_points = ret[0].size();

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					VectorFunctionType::value_type::value_type v = ret[i][qp];
					left.implementation().vector_mult(ret[i][qp], v);
				}
			}

			return ret;
		}

		template<class Op>
		inline static FunctionType apply_binary(const Number<double> &left, const TestFunction<LibMeshFunctionSpace> &right, const Op &op, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			FunctionType ret = ctx.test()[right.space_ptr()->subspace_id()]->get_phi();

			const double num = left;

			for(auto &v : ret) {
				for(auto &s : v) {
					s = op.apply(num, s);
				}
			}

			return ret;
		}

		template<class Trial, class Test>
		static Matrix bilinear_form(
			const Range &trial_range,
			const Range &test_range,
			Trial &&trial,
			Test  &&test,
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			Matrix result;
			init_tensor(result, ctx);
			
			Write<Matrix> wt(result);
			auto && dx    = ctx.dx();

			uint n_quad_points = dx.size();

			auto s = size(result);

			if(s.n_dims() == 1) {
				s.set_dims(2);
				s.set(1, 1);
			}

			const unsigned int n_test = test_range.extent();
			const unsigned int n_trial = trial_range.extent();

			for (uint i = 0; i < n_trial; i++) {
				for (uint j = 0; j < n_test; j++) {
					for (uint qp = 0; qp < n_quad_points; qp++) {
						add(result, test_range.begin() + j, trial_range.begin() + i, inner( get(trial, qp, i), get(test, qp, j) ) * dx[qp]);
					}
				}
			}

			return result;
		}

		template<class Fun, class Test>
		static Vector linear_form(
			const Range &test_range,
			Fun  &&fun,
			Test &&test,
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			Vector result;
			init_tensor(result, ctx);
			Write<Vector> wt(result);

			auto && dx    = ctx.dx();
			uint n_quad_points = dx.size();

			auto n_test = test_range.extent();
			for (uint j = 0; j < n_test; j++) {
				for (uint qp = 0; qp < n_quad_points; qp++) {
					add(result, test_range.begin() + j, 0, inner( get(fun, qp, 0), get(test, qp, j) ) * dx[qp]);
				}
			}

			return result;
		}

	};
}

#endif //UTOPIA_LIBMESH_FE_BACKEND_HPP

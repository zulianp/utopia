#ifndef UTOPIA_LIBMESH_FE_BACKEND_HPP
#define UTOPIA_LIBMESH_FE_BACKEND_HPP 

#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_AssemblyContext.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_FEBasicOverloads.hpp"

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

		template<typename T>
		inline static const T &get(const std::vector<std::vector<T> > &v, const std::size_t qp, const std::size_t i)
		{
			return v[i][qp];
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

			//Function
			// scalar fe functions
		static const FunctionType &fun(const TrialFunction<LibMeshFunctionSpace> &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return ctx.trial()[fun.space_ptr()->subspace_id()]->get_phi();
		}

		static const FunctionType &fun(const TestFunction<LibMeshFunctionSpace> &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			return ctx.test()[fun.space_ptr()->subspace_id()]->get_phi();
		}

			// vector fe functions
		static void fun_aux(
			ProductFunctionSpace<LibMeshFunctionSpace> &space,
			std::vector<std::unique_ptr<FE> > &fe_object,
			AssemblyContext<LIBMESH_TAG> &ctx,
			VectorFunctionType &ret)
		{
				// ret.resize(fe_object[0]->get_phi().size()); 

				// int n_shape_functions = 0;
				// for(auto &t : fe_object) {
				// 	n_shape_functions += t->n_shape_functions();
				// }

				// for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				// 	ret[qp].resize(n_shape_functions, zeros(space.n_subspaces()));

				// 	int offset = 0;
				// 	space.each([&offset, &fe_object, &ret, qp](const int, const LibMeshFunctionSpace &s) {
				// 		auto fe = fe_object[s.subspace_id()];

				// 		for(int j = 0; j < fe->n_shape_functions(); ++j) {
				// 			Write<Vector> w(ret[qp][offset]);

				// 			ret[qp][offset++].set(s.subspace_id(), fe->fun[qp][j]);
				// 		}
				// 	});
				// }
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
			ret.resize(fe_object[0]->get_phi().size()); 

				// int n_shape_functions = 0;
				// for(auto &t : fe_object) {
				// 	n_shape_functions += t->n_shape_functions();
				// }

				// for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				// 	ret[qp].resize(n_shape_functions, zeros(space.n_subspaces(), space[0].mesh().n_dims()));

				// 	int offset = 0;
				// 	space.each([&offset, &fe_object, &ret, &space, qp](const int, const LibMeshFunctionSpace &s) {
				// 		auto fe = fe_object[s.subspace_id()];

				// 		for(int j = 0; j < fe->n_shape_functions(); ++j, offset++) {
				// 			Write<Matrix> w(ret[qp][offset]);
				// 			Read<Vector> r(fe->grad[qp][j]);

				// 			for(int d = 0; d < s.mesh().n_dims(); ++d) {
				// 				ret[qp][offset].set(s.subspace_id(), d, fe->grad[qp][j].get(d));
				// 			}
				// 		}
				// 	});
				// }
		}

		static JacobianType grad(const TrialFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			JacobianType ret;
			grad_aux(*space_ptr, ctx.trial(), ctx, ret);
			return ret;
		}

		static JacobianType grad(const TestFunction< ProductFunctionSpace<LibMeshFunctionSpace> > &fun, AssemblyContext<LIBMESH_TAG> &ctx)
		{
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
			ret.resize(fe_object[0]->get_phi().size()); 

				// int n_shape_functions = 0;
				// for(auto &t : fe_object) {
				// 	n_shape_functions += t->n_shape_functions();
				// }

				// for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				// 	ret[qp].resize(n_shape_functions);

				// 	int offset = 0;
				// 	space.each([&offset, &fe_object, &ret, &space, qp](const int, const LibMeshFunctionSpace &s) {
				// 		auto fe = fe_object[s.subspace_id()];

				// 		for(int j = 0; j < fe->n_shape_functions(); ++j, offset++) {
				// 			Read<Vector> r(fe->grad[qp][j]);
				// 			ret[qp][offset] = fe->grad[qp][j].get(s.subspace_id()); 
				// 		}
				// 	});
				// }
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
		static void eval_curl(const Matrix &deriv, Vector &result)
		{	
				// Read<Matrix> r_x(deriv);

				// auto s = size(deriv);
				// const int dims = s.get(1);

				// switch(dims) {
				// 	case 2:
				// 	{
				// 		result = zeros(1);
				// 		Write<Vector> w_v(result);

				// 		//usually associated with the z dimension but in this case 
				// 		// we save space and save it in the x dimension
				// 		//TODO: check if there is the need to put it in the z.
				// 		result.set(0, deriv.get(1, 0) - deriv.get(0, 1));
				// 		break;
				// 	}

				// 	case 3:
				// 	{
				// 		result = zeros(3);
				// 		Write<Vector> w_v(result);

				// 		result.set(0, deriv.get(2, 1) - deriv.get(1, 2));
				// 		result.set(1, deriv.get(0, 2) - deriv.get(2, 0));
				// 		result.set(2, deriv.get(1, 0) - deriv.get(0, 1));
				// 		break;
				// 	}

				// 	default :
				// 	{
				// 		assert(false);
				// 	}
				// }
		}

		static void curl_aux(
			ProductFunctionSpace<LibMeshFunctionSpace> &space,
			std::vector<std::unique_ptr<FE> > &fe_object,
			AssemblyContext<LIBMESH_TAG> &ctx,
			CurlType &ret)
		{
				// ret.resize(fe_object[0]->get_phi().size()); 


				// JacobianType grads;

				// grad_aux(
				// 	space,
				// 	fe_object,
				// 	ctx,
				// 	grads);

				// const int n_shape_x = fe_object[0]->n_shape_functions();

				// int n_shape_functions = 0;
				// for(auto &t : fe_object) {
				// 	assert(n_shape_x == t->n_shape_functions());
				// 	n_shape_functions += t->n_shape_functions();
				// }

				// const std::size_t n_subspaces = space.n_subspaces();

				// for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				// 	ret[qp].resize(n_shape_functions);

				// 	int offset = 0;
				// 	for(int i = 0; i < n_subspaces; ++i) {
				// 		for(int j = 0; j < n_shape_x; ++j, ++offset) {

				// 			auto &ret_j  = ret[qp][offset];
				// 			auto &grad_j = grads[qp][offset];
				// 			eval_curl(grad_j, ret_j);
				// 		}
				// 	}
				// }
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
		static void interp_values(
			const Interpolate<Wrapper<Tensor, 1>, TrialFunction<LibMeshFunctionSpace> > &interp,
			Vector &element_values,
			AssemblyContext<LIBMESH_TAG> &ctx)
		{
				// auto &c   = interp.coefficient();
				// auto &f   = interp.fun();

				// auto space_ptr = f.space_ptr();

				// std::vector<int> indices;
				// space_ptr->dof_indices(ctx.current_element, indices);

				// //compute offset
				// int n_funs = ctx.trial.size();
				// for(auto &i : indices) {
				// 	i *= n_funs;
				// 	i += space_ptr->subspace_id();
				// }

				// element_values = zeros(indices.size());

				// Write<Vector> w(element_values);
				// Read<Wrapper<Tensor, 1>> r(c);

				// for(std::size_t i = 0; i < indices.size(); ++i) {
				// 	element_values.set(i, c.get(indices[i]));
				// }
		}

		template<class Tensor, class Space>
		static std::vector<Scalar> fun(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<Space> > &interp, AssemblyContext<LIBMESH_TAG> &ctx)
		{
				// auto &c  = interp.coefficient();
				// auto &f  = interp.fun();
				// auto &&g = fun(f, ctx);

				// Vector element_values;
				// interp_values(interp, element_values, ctx);

				// std::vector<Scalar> ret(g.size(), 0.);
				// for(std::size_t qp = 0; qp < g.size(); ++qp) {
				// 	for(std::size_t i = 0; i < g[qp].size(); ++i) {
				// 		ret[qp] += element_values.get(i) * g[qp][i];
				// 	}
				// }

				// return ret;
			return std::vector<Scalar>();
		}

		template<class Tensor>
		static std::vector<Vector> grad(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<LibMeshFunctionSpace> > &interp, AssemblyContext<LIBMESH_TAG> &ctx)
		{
				// auto &c   = interp.coefficient();
				// auto &f   = interp.fun();
				// auto &&g = grad(f, ctx);

				// Vector element_values;
				// interp_values(interp, element_values, ctx);

				// auto dim = size(g[0][0]).get(0);
				// std::vector<Vector> ret(g.size());
				// for(std::size_t qp = 0; qp < g.size(); ++qp) {
				// 	ret[qp] = zeros(dim);

				// 	for(std::size_t i = 0; i < g[qp].size(); ++i) {
				// 		ret[qp] += element_values.get(i) * g[qp][i];
				// 	}
				// }

				// return ret;

			return std::vector<Vector>();
		}

		template<class Tensor>
		static std::vector<Matrix> grad(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace> > > &interp, AssemblyContext<LIBMESH_TAG> &ctx)
		{
				// auto &c   = interp.coefficient();
				// auto &f   = interp.fun();
				// auto &&g = grad(f, ctx);

				// Vector element_values;
				// interp_values(interp, element_values, ctx);

				// auto s = size(g[0][0]);
				// std::vector<Matrix> ret(g.size());
				// for(std::size_t qp = 0; qp < g.size(); ++qp) {
				// 	ret[qp] = zeros(s);

				// 	for(std::size_t i = 0; i < g[qp].size(); ++i) {
				// 		ret[qp] += element_values.get(i) * g[qp][i];
				// 	}
				// }

				// return ret;
			return std::vector<Matrix>();
		}



		template<class Tensor>
		static void interp_values(
			const Interpolate<Wrapper<Tensor, 1>, TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> > &interp,
			Vector &element_values,
			AssemblyContext<LIBMESH_TAG> &ctx)
		{
				// auto &c   = interp.coefficient();
				// auto &f   = interp.fun();

				// auto space_ptr = f.space_ptr();

				// std::vector<int> indices;
				// space_ptr->subspace(0).dof_indices(ctx.current_element, indices);


				// std::vector<int> prod_indices;
				// prod_indices.reserve(indices.size() * space_ptr->n_subspaces());

				// //compute offset
				// int n_funs = ctx.trial.size();

				// space_ptr->each([&](const int sub_index, const LibMeshFunctionSpace &space) {
				// 	for(auto &i : indices) {
				// 		int p_i = i;
				// 		p_i *= n_funs;
				// 		p_i += space.subspace_id();
				// 		prod_indices.push_back(p_i);
				// 	}
				// });

				// element_values = zeros(prod_indices.size());

				// Write<Vector> w(element_values);
				// Read<Wrapper<Tensor, 1>> r(c);

				// for(std::size_t i = 0; i < prod_indices.size(); ++i) {
				// 	element_values.set(i, c.get(prod_indices[i]));
				// }
		}



			//////////////////////////////////////////////////////////////////////////////////////////

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
			JacobianType ret = grad(right, ctx);
			auto && grad_left = grad(left, ctx);

			for(std::size_t i = 0; i < ret.size(); ++i) {
				for(std::size_t qp = 0; qp < ret[i].size(); ++qp) {
					ret[i][qp] += grad_left[i][qp].transpose();
					ret[i][qp] *= scaling;
				}
			}

			return ret;
		}

			//////////////////////////////////////////////////////////////////////////////////////////////////////

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

		// template<class Space>
		inline static auto multiply(
			const LMDenseMatrix &left,
			const Gradient<TrialFunction<ProductFunctionSpace<LibMeshFunctionSpace>> > &right, 
			AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type 
		{
			typename remove_ref_and_const<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			Read<LMDenseMatrix> r_l(left);

			Size s_l = size(left);

			for(auto &v : ret) {
				for(auto &s : v) {
					auto s_copy = s;

					for(uint k = 0; k < LIBMESH_DIM; ++k) {
						for(uint i = 0; i < s_l.get(0); ++i) {
							s(i, k) = left.get(i, 0) * s_copy(0, k);

							for(uint j = 1; j < s_l.get(1); ++j) {
								s(i, k) += left.get(i, j) * s_copy(j, k);
							}	
						}
					}
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

			const std::size_t n_functions = right.size();
			const std::size_t n_quad_points = right[0].size();

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
			const std::size_t n_functions = right.size();

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					ret[i][qp] = left[qp] * ret[i][qp];
				}
			}

			return ret;
		}

		template<class Space>
		inline static auto multiply(const std::vector<Matrix> &left, const TrialFunction<Space> &right, AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(fun(right, ctx))>::type
		{
			typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			const std::size_t n_quad_points = left.size();
			const std::size_t n_functions = right.size();

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					ret[i][qp] = left[qp] * ret[i][qp];
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
		inline static auto multiply(const Left &left, const TestFunction<Space> &right,  AssemblyContext<LIBMESH_TAG> &ctx) -> typename remove_ref_and_const<decltype(fun(right, ctx))>::type
		{
			typename remove_ref_and_const<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			const std::size_t n_functions = right.size();
			const std::size_t n_quad_points = right[0].size();

			for(std::size_t qp = 0; qp < n_quad_points; ++qp) {
				for(std::size_t i = 0; i < n_functions; ++i) {
					ret[i][qp] = left * ret[i][qp];
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
			Trial &&trial,
			Test  &&test,
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			Matrix result;
			ctx.init_tensor(result, true);
			Write<Matrix> wt(result);

			auto && dx    = ctx.dx();

			uint n_quad_points = dx.size();

			auto s = size(result);

			if(s.n_dims() == 1) {
				s.set_dims(2);
				s.set(1, 1);
			}

			for (uint i = 0; i < s.get(1); i++) {
				for (uint j = 0; j < s.get(0); j++) {
					for (uint qp = 0; qp < n_quad_points; qp++) {
						add(result, j, i, inner( get(trial, qp, i), get(test, qp, j) ) * dx[qp]);
					}
				}
			}

			return result;
		}

		template<class Fun, class Test>
		static Vector linear_form(
			Fun  &&fun,
			Test &&test,
			AssemblyContext<LIBMESH_TAG> &ctx)
		{	
			Vector result;
			ctx.init_tensor(result, true);
			Write<Vector> wt(result);

			auto && dx    = ctx.dx();
			uint n_quad_points = dx.size();

			auto s = size(result);
			for (uint j = 0; j < s.get(0); j++) {
				for (uint qp = 0; qp < n_quad_points; qp++) {
					add(result, j, 0, inner( get(fun, qp, 0), get(test, qp, j) ) * dx[qp]);
				}
			}

			return result;
		}
	};
}

#endif //UTOPIA_LIBMESH_FE_BACKEND_HPP

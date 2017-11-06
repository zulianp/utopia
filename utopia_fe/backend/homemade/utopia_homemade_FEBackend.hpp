#ifndef UTOPIA_HOMEADE_FE_BACKEND_HPP
#define UTOPIA_HOMEADE_FE_BACKEND_HPP 

#include "utopia_Base.hpp"
#include "utopia_homemade_FEForwardDeclarations.hpp"

#include "utopia_FEBackend.hpp"
#include "utopia_fe_core.hpp"

#include "utopia_homemade_FunctionSpace.hpp"
#include <vector>

namespace utopia {


	inline static double inner(const double left, const double right)
	{
		return left * right;
	}

	template<>
	class FEBackend<HOMEMADE> {
	public:
		template<typename T>
		inline static const T &get(const std::vector<std::vector<T> > &v, const std::size_t qp, const std::size_t i)
		{
			return v[qp][i];
		}

		template<typename T, int Order>
		inline static const T get(const ConstantCoefficient<T, Order> &c, const std::size_t qp, const std::size_t i)
		{
			return c[i];
		}

		template<typename T>
		inline static void add(ElementMatrix &mat, const int i, const int j, const T value)
		{
			mat.add(i, j, value);
		}

		template<typename T>
		inline static void add(ElementVector &vec, const int i, const int j, const T value)
		{
			assert(j == 0);
			vec.add(i, value);
		}

		//////////////////////////////////////////////////////////////////////////////////////////

		//Function
		// scalar fe functions
		static HMFun &fun(const TrialFunction<HMFESpace> &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			return ctx.trial[fun.space_ptr()->subspace_id()]->fun;
		}

		static HMFun &fun(const TestFunction<HMFESpace> &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			return ctx.test[fun.space_ptr()->subspace_id()]->fun;
		}

		// vector fe functions
		static void fun_aux(
			ProductFunctionSpace<HMFESpace> &space,
			std::vector<std::shared_ptr<FE> > &fe_object,
			AssemblyContext<HOMEMADE> &ctx,
			HMVectorFEFun &ret)
		{
			ret.resize(fe_object[0]->fun.size()); 

			int n_shape_functions = 0;
			for(auto &t : fe_object) {
				n_shape_functions += t->n_shape_functions();
			}

			for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				ret[qp].resize(n_shape_functions, zeros(space.n_subspaces()));

				int offset = 0;
				space.each([&offset, &fe_object, &ret, qp](const int, const HMFESpace &s) {
					auto fe = fe_object[s.subspace_id()];
					
					for(int j = 0; j < fe->n_shape_functions(); ++j) {
						Write<Vectord> w(ret[qp][offset]);

						ret[qp][offset++].set(s.subspace_id(), fe->fun[qp][j]);
					}
				});
			}
		}

		static HMVectorFEFun fun(const TestFunction< ProductFunctionSpace<HMFESpace> > &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			HMVectorFEFun ret;

			fun_aux(*space_ptr, ctx.test, ctx, ret);
			return ret;
		}

		static HMVectorFEFun fun(const TrialFunction< ProductFunctionSpace<HMFESpace> > &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			HMVectorFEFun ret;
			fun_aux(*space_ptr, ctx.trial, ctx, ret);
			return ret;
		}

		//////////////////////////////////////////////////////////////////////////////////////////

		//Gradient
		// scalar fe functions
		static HMDerivative &grad(const TrialFunction<HMFESpace> &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			return ctx.trial[fun.space_ptr()->subspace_id()]->grad;
		}

		static HMDerivative &grad(const TestFunction<HMFESpace> &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			return ctx.test[fun.space_ptr()->subspace_id()]->grad;
		}


		// vector fe functions
		static void grad_aux(
			ProductFunctionSpace<HMFESpace> &space,
			std::vector<std::shared_ptr<FE> > &fe_object,
			AssemblyContext<HOMEMADE> &ctx,
			HMVectorFEDerivative &ret)
		{
			ret.resize(fe_object[0]->fun.size()); 

			int n_shape_functions = 0;
			for(auto &t : fe_object) {
				n_shape_functions += t->n_shape_functions();
			}

			for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				ret[qp].resize(n_shape_functions, zeros(space.n_subspaces(), space[0].mesh().n_dims()));

				int offset = 0;
				space.each([&offset, &fe_object, &ret, &space, qp](const int, const HMFESpace &s) {
					auto fe = fe_object[s.subspace_id()];
					
					for(int j = 0; j < fe->n_shape_functions(); ++j, offset++) {
						Write<Matrixd> w(ret[qp][offset]);
						Read<Vectord> r(fe->grad[qp][j]);

						for(int d = 0; d < s.mesh().n_dims(); ++d) {
							ret[qp][offset].set(s.subspace_id(), d, fe->grad[qp][j].get(d));
						}
					}
				});
			}
		}

		static HMVectorFEDerivative grad(const TrialFunction< ProductFunctionSpace<HMFESpace> > &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			HMVectorFEDerivative ret;
			grad_aux(*space_ptr, ctx.trial, ctx, ret);
			return ret;
		}

		static HMVectorFEDerivative grad(const TestFunction< ProductFunctionSpace<HMFESpace> > &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			HMVectorFEDerivative ret;
			grad_aux(*space_ptr, ctx.test, ctx, ret);
			return ret;
		}


		//////////////////////////////////////////////////////////////////////////////////////////

		//Divergence
		// vector fe functions
		static void div_aux(
			ProductFunctionSpace<HMFESpace> &space,
			std::vector<std::shared_ptr<FE> > &fe_object,
			AssemblyContext<HOMEMADE> &ctx,
			HMFun &ret)
		{
			ret.resize(fe_object[0]->fun.size()); 

			int n_shape_functions = 0;
			for(auto &t : fe_object) {
				n_shape_functions += t->n_shape_functions();
			}

			for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				ret[qp].resize(n_shape_functions);

				int offset = 0;
				space.each([&offset, &fe_object, &ret, &space, qp](const int, const HMFESpace &s) {
					auto fe = fe_object[s.subspace_id()];
					
					for(int j = 0; j < fe->n_shape_functions(); ++j, offset++) {
						Read<Vectord> r(fe->grad[qp][j]);
						ret[qp][offset] = fe->grad[qp][j].get(s.subspace_id()); 
					}
				});
			}
		}

		static HMFun div(const TrialFunction< ProductFunctionSpace<HMFESpace> > &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			HMFun ret;
			div_aux(*space_ptr, ctx.trial, ctx, ret);
			return ret;
		}

		static HMFun div(const TestFunction< ProductFunctionSpace<HMFESpace> > &fun, AssemblyContext<HOMEMADE> &ctx)
		{
			auto space_ptr = fun.space_ptr();
			HMFun ret;
			div_aux(*space_ptr, ctx.test, ctx, ret);
			return ret;
		}

		//////////////////////////////////////////////////////////////////////////////////////////
		//Interpolate
		template<class Tensor>
		static void interp_values(
			const Interpolate<Wrapper<Tensor, 1>, TrialFunction<HMFESpace> > &interp,
			Vectord &element_values,
			AssemblyContext<HOMEMADE> &ctx)
		{
			auto &c   = interp.coefficient();
			auto &f   = interp.fun();

			auto space_ptr = f.space_ptr();

			std::vector<int> indices;
			space_ptr->dof_indices(ctx.current_element, indices);

			//compute offset
			int n_funs = ctx.trial.size();
			for(auto &i : indices) {
			 	 i *= n_funs;
			 	 i += space_ptr->subspace_id();
			}
			
			element_values = zeros(indices.size());

			Write<Vectord> w(element_values);
			Read<Wrapper<Tensor, 1>> r(c);

			for(std::size_t i = 0; i < indices.size(); ++i) {
				element_values.set(i, c.get(indices[i]));
			}
		}

		template<class Tensor, class Space>
		static std::vector<double> fun(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<Space> > &interp, AssemblyContext<HOMEMADE> &ctx)
		{
			auto &c  = interp.coefficient();
			auto &f  = interp.fun();
			auto &&g = fun(f, ctx);

			Vectord element_values;
			interp_values(interp, element_values, ctx);
			
			std::vector<double> ret(g.size(), 0.);
			for(std::size_t qp = 0; qp < g.size(); ++qp) {
				for(std::size_t i = 0; i < g[qp].size(); ++i) {
					ret[qp] += element_values.get(i) * g[qp][i];
				}
			}

			return ret;
		}

		template<class Tensor>
		static std::vector<Vectord> grad(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<HMFESpace> > &interp, AssemblyContext<HOMEMADE> &ctx)
		{
			auto &c   = interp.coefficient();
			auto &f   = interp.fun();
			auto &&g = grad(f, ctx);

			Vectord element_values;
			interp_values(interp, element_values, ctx);
			
			auto dim = size(g[0][0]).get(0);
			std::vector<Vectord> ret(g.size());
			for(std::size_t qp = 0; qp < g.size(); ++qp) {
				ret[qp] = zeros(dim);

				for(std::size_t i = 0; i < g[qp].size(); ++i) {
					ret[qp] += element_values.get(i) * g[qp][i];
				}
			}

			return ret;
		}

		template<class Tensor>
		static std::vector<Matrixd> grad(const Interpolate<Wrapper<Tensor, 1>, TrialFunction<ProductFunctionSpace<HMFESpace> > > &interp, AssemblyContext<HOMEMADE> &ctx)
		{
			auto &c   = interp.coefficient();
			auto &f   = interp.fun();
			auto &&g = grad(f, ctx);

			Vectord element_values;
			interp_values(interp, element_values, ctx);
			
			auto s = size(g[0][0]);
			std::vector<Matrixd> ret(g.size());
			for(std::size_t qp = 0; qp < g.size(); ++qp) {
				ret[qp] = zeros(s);

				for(std::size_t i = 0; i < g[qp].size(); ++i) {
					ret[qp] += element_values.get(i) * g[qp][i];
				}
			}

			return ret;
		}



		template<class Tensor>
		static void interp_values(
			const Interpolate<Wrapper<Tensor, 1>, TrialFunction<ProductFunctionSpace<HMFESpace>> > &interp,
			Vectord &element_values,
			AssemblyContext<HOMEMADE> &ctx)
		{
			auto &c   = interp.coefficient();
			auto &f   = interp.fun();

			auto space_ptr = f.space_ptr();

			std::vector<int> indices;
			space_ptr->subspace(0).dof_indices(ctx.current_element, indices);


			std::vector<int> prod_indices;
			prod_indices.reserve(indices.size() * space_ptr->n_subspaces());

			//compute offset
			int n_funs = ctx.trial.size();

			space_ptr->each([&](const int sub_index, const HMFESpace &space) {
				for(auto &i : indices) {
					 int p_i = i;
				 	 p_i *= n_funs;
				 	 p_i += space.subspace_id();
				 	 prod_indices.push_back(p_i);
				}
			});
			
			element_values = zeros(prod_indices.size());

			Write<Vectord> w(element_values);
			Read<Wrapper<Tensor, 1>> r(c);

			for(std::size_t i = 0; i < prod_indices.size(); ++i) {
				element_values.set(i, c.get(prod_indices[i]));
			}
		}



		//////////////////////////////////////////////////////////////////////////////////////////

		template<class Space, class Op>
		inline static auto apply_binary(const Number<double> &left, const TrialFunction<Space> &right, const Op &op, AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(fun(right, ctx))>::type 
		{
			typename std::remove_reference<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			const double num = left;
			for(auto &v : ret) {
				for(auto &s : v) {
					s = op.apply(num, s);
				}
			}

			return ret;
		}

		template<class Space, class Op>
		inline static auto apply_binary(const Number<double> &left, const TestFunction<Space> &right, const Op &op, AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(fun(right, ctx))>::type 
		{
			typename std::remove_reference<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

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
		inline static HMVectorFEDerivative grad_t_plus_grad(const double scaling, const Left &left, const Right &right, AssemblyContext<HOMEMADE> &ctx)
		{
			HMVectorFEDerivative ret = grad(right, ctx);
			auto && grad_left = grad(left, ctx);

			for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				for(std::size_t i = 0; i < ret[qp].size(); ++i) {
					ret[qp][i] += transpose(grad_left[qp][i]);
					ret[qp][i] *= scaling;
				}
			}

			return ret;
		}


		//////////////////////////////////////////////////////////////////////////////////////////////////////

		template<class Left, class Space>
		inline static auto multiply(
			const Left &left,
			const Gradient<TrialFunction<Space> > &right, 
			AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(grad(right.expr(), ctx))>::type 
		{
			// HMDerivative ret = ctx.trial[right.expr().space_ptr()->subspace_id()]->grad;
			typename std::remove_reference<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			for(auto &v : ret) {
				for(auto &s : v) {
					s = left * s;
				}
			}

			return ret;
		}


		template<class Space>
		inline static auto multiply(
			const std::vector<double> &left,
			const Gradient<TrialFunction<Space> > &right, 
			AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(grad(right.expr(), ctx))>::type 
		{
			typename std::remove_reference<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			for(std::size_t i = 0; i < ret.size(); ++i) {
				auto &v = ret[i];

				for(auto &s : v) {
					s = left[i] * s;
				}
			}

			return ret;
		}





		template<class Left, class Space>
		inline static auto multiply(
			const Left &left,
			const Gradient<TestFunction<Space> > &right, 
			AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(grad(right.expr(), ctx))>::type 
		{
			// HMDerivative ret = ctx.trial[right.expr().space_ptr()->subspace_id()]->grad;
			typename std::remove_reference<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			for(auto &v : ret) {
				for(auto &s : v) {
					s = left * s;
				}
			}

			return ret;
		}

		template<class Space>
		inline static auto multiply(const std::vector<double> &left, const TrialFunction<Space> &right,  AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(fun(right, ctx))>::type
		{
			typename std::remove_reference<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				auto &v = ret[qp];
				for(auto &s : v) {
					s = left[qp] * s;
				}
			}

			return ret;
		}


		template<class Space>
		inline static auto multiply(const std::vector<Matrixd> &left, const Gradient<TrialFunction<Space> > &right,  AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(grad(right.expr(), ctx))>::type
		{
			typename std::remove_reference<decltype(grad(right.expr(), ctx))>::type ret = grad(right.expr(), ctx);

			for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				auto &v = ret[qp];
				for(auto &s : v) {
					s = left[qp] * s;
				}
			}

			return ret;
		}

		template<class Space>
		inline static auto multiply(const std::vector<Matrixd> &left, const TrialFunction<Space> &right,  AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(fun(right, ctx))>::type
		{
			typename std::remove_reference<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			for(std::size_t qp = 0; qp < ret.size(); ++qp) {
				auto &v = ret[qp];
				for(auto &s : v) {
					s = left[qp] * s;
				}
			}

			return ret;
		}

		template<class Left, class Space>
		inline static auto multiply(const Left &left, const TrialFunction<Space> &right,  AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(fun(right, ctx))>::type
		{
			typename std::remove_reference<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			for(auto &v : ret) {
				for(auto &s : v) {
					s = left * s;
				}
			}

			return ret;
		}

		template<class Left, class Space>
		inline static auto multiply(const Left &left, const TestFunction<Space> &right,  AssemblyContext<HOMEMADE> &ctx) -> typename std::remove_reference<decltype(fun(right, ctx))>::type
		{
			typename std::remove_reference<decltype(fun(right, ctx))>::type ret = fun(right, ctx);

			for(auto &v : ret) {
				for(auto &s : v) {
					s = left * s;
				}
			}

			return ret;
		}

		template<class Op>
		inline static HMFun apply_binary(const Number<double> &left, const TestFunction<HMFESpace> &right, const Op &op, AssemblyContext<HOMEMADE> &ctx)
		{
			HMFun ret = ctx.test[right.space_ptr()->subspace_id()]->fun;
			
			const double num = left;
			
			for(auto &v : ret) {
				for(auto &s : v) {
					s = op.apply(num, s);
				}
			}

			return ret;
		}

		template<class Trial, class Test>
		static Matrixd bilinear_form(
		 	Trial &&trial,
			Test  &&test,
			AssemblyContext<HOMEMADE> &ctx)
		{	
			Matrixd result;
			ctx.init_tensor(result, true);
			Write<Matrixd> wt(result);

			auto && dx    = ctx.dx();


			uint n_quad_points = dx.size();

			auto s = size(result);

			if(s.n_dims() == 1) {
				s.set_dims(2);
				s.set(1, 1);
			}


			for (uint qp = 0; qp < n_quad_points; qp++) {
				for (uint i = 0; i < s.get(1); i++) {
					for (uint j = 0; j < s.get(0); j++) {
						add(result, j, i, inner( get(trial, qp, i), get(test, qp, j) ) * dx[qp]);
					}
				}
			}

			return result;
		}

		template<class Fun, class Test>
		static Vectord linear_form(
		 	Fun  &&fun,
			Test &&test,
			AssemblyContext<HOMEMADE> &ctx)
		{	
			Vectord result;
			ctx.init_tensor(result, true);
			Write<Vectord> wt(result);

			auto && dx    = ctx.dx();
			uint n_quad_points = dx.size();

			auto s = size(result);
			for (uint qp = 0; qp < n_quad_points; qp++) {
				for (uint j = 0; j < s.get(0); j++) {
					add(result, j, 0, inner( get(fun, qp, 0), get(test, qp, j) ) * dx[qp]);
				}
			}

			return result;
		}

	};

}

#endif //UTOPIA_HOMEADE_FE_BACKEND_HPP

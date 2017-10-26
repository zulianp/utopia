#ifndef UTOPIA_FE_HOMEMADE_HPP
#define UTOPIA_FE_HOMEMADE_HPP 


#include "utopia_Traits.hpp"
#include "utopia_Base.hpp"
#include "utopia_fe_core.hpp"
#include "utopia_intersector.hpp"
#include "utopia_FEBackend.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FormEval.hpp"
#include "utopia_FEEval.hpp"

#include <memory>

namespace utopia {

	typedef utopia::Matrixd ElementMatrix;
	typedef utopia::Vectord ElementVector;
	typedef std::vector< std::vector<utopia::ElementVector> > HMDerivative;
	typedef std::vector< std::vector<utopia::ElementMatrix> > HMJacobian;
	typedef std::vector< std::vector<double> > HMFun;

	typedef std::vector< std::vector<utopia::ElementVector> > HMVectorFEFun;
	typedef std::vector< std::vector<utopia::ElementMatrix> > HMVectorFEDerivative;

	typedef std::vector<double> HMDx;


	static const int EMPTY_BACKEND_FLAG = -3;

	class Empty {};

	template<>
	class Traits<Empty> {
	public:
		static const int Backend = EMPTY_BACKEND_FLAG;
		static const int Order = 1;
		typedef double Scalar;
		typedef int SizeType;
		typedef utopia::Empty Implementation;

	};

	class EmptyFunctionSpace : public FunctionSpace<EmptyFunctionSpace> {
	public:
	};

	template<>
	class FormTraits<EmptyFunctionSpace> : public Traits<Empty> {};

	class Mesh {
	public:
		typedef Intersector::Mesh MeshImpl;

		std::shared_ptr<MeshImpl> impl;

		void make_example_mesh()
		{
			impl = std::make_shared<MeshImpl>();

			el_ptr.resize(2);
			el_ptr[0] = 0;
			el_ptr[1] = 3;

			el_index.resize(3);
			el_index[0] = 0;
			el_index[1] = 1;
			el_index[2] = 2;

			el_type.resize(1);
			el_type[0] = Intersector::ELEMENT_TYPE_TRIANGLE;

			points.resize(3 * 2);
			points[0] = 0.0;
			points[1] = 0.0;

			points[2] = 1.0;
			points[3] = 0.0;

			points[4] = 0.0;
			points[5] = 1.0;

			impl->n_elements = 1;
			impl->n_nodes  = 3;
			impl->n_dims   = 2;
			impl->el_ptr   = &el_ptr[0];
			impl->el_index = &el_index[0];
			impl->el_type  = &el_type[0];
			impl->meta     = nullptr;// &meta[0];
			impl->points   = &points[0];
		}

		Mesh()
		{
			make_example_mesh();
		}

	private:
		//memory
		std::vector<int> el_ptr;
		std::vector<int> el_index;
		std::vector<int> el_type;
		std::vector<int> meta;
		std::vector<double> points;
	};

	class HMFESpace : public FunctionSpace<HMFESpace> {
	public:
		Mesh mesh;
	};

	template<>
	class Traits<HMFESpace> {
	public:
		static const int Backend = HOMEMADE;
		static const int Order = 1;
		typedef double Scalar;
		typedef utopia::HMFESpace Implementation;
		typedef utopia::HMDerivative GradientType;
		typedef utopia::HMJacobian JacobianType;
	};

	class Quadrature {
	public:
		double * points() { return nullptr; }
		int n_points() { return 0; }
	};

	class FE {
	public:
		void init(const int current_element, Mesh &mesh, int quadrature_order)
		{

			// fe_backend.make_fe_object_from_quad_2(
			// 	*mesh.impl, current_element, 3, quad_points, quad_weights, &fe_);
			fe_backend.make_fe_object(*mesh.impl, current_element, quadrature_order, &fe_);

			grad.resize(fe_.n_quad_points);
			fun.resize(grad.size());
			dx.resize(grad.size());

			for(int q = 0; q < fe_.n_quad_points; ++q) {
				grad[q].resize(fe_.n_shape_functions);
				fun[q].resize(fe_.n_shape_functions);

				dx[q] = fe_.dx[q];

				for(int i = 0; i < fe_.n_shape_functions; ++i) {
					grad[q][i] = zeros(fe_.n_dims);

					Write<ElementVector> w(grad[q][i]);
					fun[q][i] = fe_.fun[i][q];
					
					for(int k = 0; k < fe_.n_dims; ++k) {
						grad[q][i].set(k, fe_.grad[i][q * fe_.n_dims + k]);
					}
				}
			}
		}

		HMDerivative grad;
		HMFun fun;
		HMDx dx;

		inline int n_shape_functions() const
		{
			return fe_.n_shape_functions;
		}

	private:

		Intersector fe_backend;
		FEObject fe_;
		
	};

	template<>
	class AssemblyContext<HOMEMADE> {
	public:

		std::vector< std::shared_ptr<FE> > trial;
		std::vector< std::shared_ptr<FE> > test;

		long quadrature_order;
		long current_element;

		template<class Expr>
		void init_bilinear(const Expr &expr) 
		{		
			//clean up previous context
			trial.clear();
			test.clear();

			std::shared_ptr<FE> trial_fe;
			std::shared_ptr<FE> test_fe;

			auto trial_space_ptr = trial_space<HMFESpace>(expr);
			auto test_space_ptr  = test_space<HMFESpace>(expr);

			if(trial_space_ptr) {
				trial_fe = std::make_shared<FE>();
				trial_fe->init(current_element, trial_space_ptr->mesh, quadrature_order);
				trial.push_back(trial_fe);
			} else {
				//product space init
				auto prod_trial_space_ptr = trial_space<ProductFunctionSpace<HMFESpace>>(expr);
				assert(prod_trial_space_ptr);

				trial_fe = std::make_shared<FE>();
				trial_fe->init(current_element, prod_trial_space_ptr->subspace(0).mesh, quadrature_order);

				//just copy it three times for now
				for(std::size_t i = 0; i < prod_trial_space_ptr->n_subspaces(); ++i) {
					trial.push_back(trial_fe);
				}
			}

			if(test_space_ptr) {
				if(trial_space_ptr != test_space_ptr) {
					test_fe = std::make_shared<FE>();
					test_fe->init(current_element, test_space_ptr->mesh, quadrature_order);
					test.push_back(test_fe);
				} else if(trial_space_ptr) {
					test_fe = trial_fe;
				}
			} else {
				//product space init
				auto prod_test_space_ptr = test_space<ProductFunctionSpace<HMFESpace>>(expr);
				assert(prod_test_space_ptr);

				test_fe = std::make_shared<FE>();
				test_fe->init(current_element, prod_test_space_ptr->subspace(0).mesh, quadrature_order);

				//just copy it three times for now
				for(std::size_t i = 0; i < prod_test_space_ptr->n_subspaces(); ++i) {
					test.push_back(test_fe);
				}
			}
		}

		template<class Expr>
		void init_linear(const Expr &expr) 
		{		
			test.clear();
			
			auto test_space_ptr  = test_space<HMFESpace>(expr);
			std::shared_ptr<FE> test_fe;

			if(test_space_ptr) {
				test_fe = std::make_shared<FE>();
				test_fe->init(current_element, test_space_ptr->mesh, quadrature_order);
				test.push_back(test_fe);
			
			} else {
				//product space init
				auto prod_test_space_ptr = test_space<ProductFunctionSpace<HMFESpace>>(expr);
				assert(prod_test_space_ptr);

				test_fe = std::make_shared<FE>();
				test_fe->init(current_element, prod_test_space_ptr->subspace(0).mesh, quadrature_order);

				//just copy it three times for now
				for(std::size_t i = 0; i < prod_test_space_ptr->n_subspaces(); ++i) {
					test.push_back(test_fe);
				}
			}
		}

		void init_tensor(ElementVector &v, const bool reset) {
			auto s = size(v);

			int n_shape_functions = 0;
			for(auto &t : test) {
				n_shape_functions += t->n_shape_functions();
			}

			if(reset || s.get(0) != n_shape_functions) {
				v = zeros(n_shape_functions);
			}
		}

		void init_tensor(ElementMatrix &v, const bool reset) {
			auto s = size(v);

			int n_trial_functions = 0;
			for(auto &t : trial) {
				n_trial_functions += t->n_shape_functions();
			}

			int n_test_functions = 0;
			for(auto &t : test) {
				n_test_functions += t->n_shape_functions();
			}

			if(reset || s.get(0) != n_test_functions || s.get(1) != n_trial_functions) {
				v = zeros(n_test_functions, n_trial_functions);
			}
		}

		const HMDx &dx() const
		{
			return test[0]->dx;
		}

		AssemblyContext()
		: current_element(0), quadrature_order(2)
		{}
	};

	template<>
	class FEBackend<HOMEMADE> {
	public:

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
				ret[qp].resize(n_shape_functions, zeros(space.n_subspaces(), space[0].mesh.impl->n_dims));

				int offset = 0;
				space.each([&offset, &fe_object, &ret, &space, qp](const int, const HMFESpace &s) {
					auto fe = fe_object[s.subspace_id()];
					
					for(int j = 0; j < fe->n_shape_functions(); ++j, offset++) {
						Write<Matrixd> w(ret[qp][offset]);
						Read<Vectord> r(fe->grad[qp][j]);

						for(int d = 0; d < s.mesh.impl->n_dims; ++d) {
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

	};

	inline static double inner(const double left, const double right)
	{
		return left * right;
	}

	template<class Form>
	class FormEval<Form, HOMEMADE> {
	public:

		typedef utopia::Traits<HMFESpace> Traits;
		FormEval() { }

		template<class Expr, class Tensor>
		static void apply(
			const Integral<Expr> &expr, 
			Tensor &mat, 
			AssemblyContext<HOMEMADE> &ctx)
		{
			apply(expr.expr(), mat, ctx);
		}

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

		template<class Left, class Right, class Tensor>
		static void apply(
			const Reduce<Binary<Left, Right, EMultiplies>, Plus> &expr, 
			Tensor &result, 
			AssemblyContext<HOMEMADE> &ctx)
		{	
			Write<Tensor> wt(result);

			auto && left  = FEEval<Left,  Traits, HOMEMADE>::apply(expr.expr().left(),  ctx);
			auto && right = FEEval<Right, Traits, HOMEMADE>::apply(expr.expr().right(), ctx);
			auto && dx    = ctx.dx();

			const bool left_is_test = is_test(expr.expr().left());
			assert( left_is_test != is_test(expr.expr().right()) );

			uint n_quad_points = dx.size();

			auto s = size(result);
			
			if(s.n_dims() == 1) {
				s.set_dims(2);
				s.set(1, 1);
			}
			
			if(left_is_test) {
				for (uint qp = 0; qp < n_quad_points; qp++) {
					for (uint i = 0; i < s.get(0); i++) {
						for (uint j = 0; j < s.get(1); j++) {
							add(result, i, j, inner( get(left, qp, i), get(right, qp, j) ) * dx[qp]);
						}
					}
				}

			} else {
				for (uint qp = 0; qp < n_quad_points; qp++) {
					for (uint i = 0; i < s.get(1); i++) {
						for (uint j = 0; j < s.get(0); j++) {
							add(result, j, i, inner( get(left, qp, i), get(right, qp, j) ) * dx[qp]);
						}
					}
				}
			}
		}
	};
}

#endif //UTOPIA_FE_HOMEMADE_HPP

#ifndef UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP
#define UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP 

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_Types.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_FunctionalTraits.hpp"
#include "utopia_libmesh_Utils.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_Traverse.hpp"
#include "utopia_ProductFunctionSpace.hpp"

#include "utopia_fe_core.hpp"

#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/fe_interface.h"

namespace utopia {

	class LibMeshAssemblyContext {
	public:
		typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
		typedef TraitsT::FE FE; 
		typedef TraitsT::Matrix Matrix;
		typedef TraitsT::Vector Vector;
		typedef TraitsT::DXType DXType;

		inline std::vector< std::unique_ptr<FE> > &fe()		
		{
			return fe_;
		}


		inline const std::vector< std::unique_ptr<FE> > &fe() const
		{
			return fe_;
		}

		inline std::vector< std::unique_ptr<FE> > &test()		
		{
			return fe_;
		}

		inline std::vector< std::unique_ptr<FE> > &trial()
		{
			return fe_;
		}

		inline const std::vector< std::unique_ptr<FE> > &test()	const
		{
			return fe_;
		}

		inline const std::vector< std::unique_ptr<FE> > &trial() const
		{
			return fe_;
		}

		inline std::shared_ptr<libMesh::QBase> quad_test() const
		{
			return quad_test_;
		}

		inline std::shared_ptr<libMesh::QBase> quad_trial() const
		{
			if(!quad_trial_) return quad_test_;
			return quad_trial_;
		}

		template<class Expr>
		void init_bilinear(const Expr &expr) 
		{		
			static_assert( (IsSubTree<TrialFunction<utopia::Any>, Expr>::value), "could not find trial function" );
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function"  );
			init_fe_from(expr);
			// std::cout << "init_bilinear: quadrature_order: " << quadrature_order_ << std::endl;
		}

		template<class Expr>
		void init_linear(const Expr &expr) 
		{	
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function"  );	
			init_fe_from(expr);
			// std::cout << "init_linear: quadrature_order: " << quadrature_order_ << std::endl;
		}

		template<class Expr>
		static std::shared_ptr<LibMeshFunctionSpace> find_any_space(const Expr &expr)
		{
			std::shared_ptr<LibMeshFunctionSpace> space_ptr = test_space<LibMeshFunctionSpace>(expr);
			if(!space_ptr) {
				auto prod_space_ptr = test_space<ProductFunctionSpace<LibMeshFunctionSpace>>(expr);
				space_ptr = prod_space_ptr->subspace_ptr(0);
			}

			return space_ptr;
		}

		template<class Expr>
		void init_offsets(const Expr &expr) 
		{
			auto space_ptr = find_any_space(expr);

			assert((static_cast<bool>(space_ptr)));

			const std::size_t n_vars = space_ptr->equation_system().n_vars();
			offset.resize(n_vars + 1);
			offset[0] = 0;

			const libMesh::Elem * elem = space_ptr->mesh().elem(current_element_);

			const libMesh::ElemType type = elem->type();
			const unsigned int sys_num   = space_ptr->dof_map().sys_number();
			const unsigned int dim       = elem->dim();

			for(std::size_t i = 0; i < n_vars; ++i) {
				const auto & var = space_ptr->dof_map().variable(i);
				libMesh::FEType fe_type = var.type();
				offset[i+1] = offset[i] + libMesh::FEInterface::n_shape_functions(dim, fe_type, type);
			}
		}

		template<class Expr>		
		void init_fe_flags(const Expr &expr)
		{
			FEInitializer fe_init(*this);
			fe_init.apply(expr);
		}

		template<class Expr>
		void init_fe_from(const Expr &expr)
		{
			auto space_ptr = find_any_space(expr);
			space_ptr->initialize();
			quadrature_order_ = functional_order(expr, *this);
			const int dim = space_ptr->mesh().mesh_dimension();
			const libMesh::Elem * elem = space_ptr->mesh().elem(current_element_);

			if(is_quad(elem->type())) {
				const int temp = quadrature_order_/2;
				quadrature_order_ = (temp + 1) * 2;
			} else if(is_hex(elem->type())) {
				const int temp = quadrature_order_/2;
				quadrature_order_ = (temp + 2) * 2;
			}

			init_offsets(expr);
			

			set_up_quadrature(dim, quadrature_order_, elem);
			block_id_ = elem->subdomain_id();

			const auto &eq_sys = space_ptr->equation_system();
			const std::size_t n_vars = eq_sys.n_vars();
			fe_.resize(n_vars);

			for(std::size_t i = 0; i < n_vars; ++i) {
				auto fe = libMesh::FEBase::build(dim, eq_sys.get_dof_map().variable_type(i));
				fe->attach_quadrature_rule(quad_test().get());
				fe_[i] = std::move(fe);
						
			}

			init_fe_flags(expr);

			for(std::size_t i = 0; i < n_vars; ++i) {
				fe_[i]->reinit(elem);
			}
		}

		void init_tensor(Vector &v, const bool reset);
		void init_tensor(Matrix &v, const bool reset);

		const DXType &dx() const
		{
			return test()[0]->get_JxW();
		}

		inline int block_id() const
		{
			return block_id_;
		}

		inline void set_current_element(const long current_element)
		{
			current_element_ = current_element;
		}

		inline long current_element() const
		{
			return current_element_;
		}

		void surface_integral_begin()
		{
			is_surface_ = true;
		}

		void surface_integral_end()
		{
			is_surface_ = false;
		}

		LibMeshAssemblyContext()
		: current_element_(0), quadrature_order_(2), block_id_(0), reset_quadrature_(true), is_surface_(false)
		{}

		//x basis function
		std::vector<int> offset;
	private:
		long quadrature_order_;
		long current_element_;
		int block_id_;
		bool reset_quadrature_;
		bool is_surface_;
		
		std::shared_ptr<libMesh::QBase> quad_trial_;
		std::shared_ptr<libMesh::QBase> quad_test_;
		std::vector< std::unique_ptr<FE> > fe_;

		inline std::size_t n_shape_functions() const
		{
			std::size_t ret = 0;
			for(auto &f_ptr : fe()) {
				if(f_ptr) ret += f_ptr->n_shape_functions();
			}

			return ret;
		}

		void set_up_quadrature(const int dim, const int quadrature_order, const libMesh::Elem * elem)
		{
			if(reset_quadrature_) {
				quad_test_  = std::make_shared<libMesh::QGauss>(dim, libMesh::Order(quadrature_order));
				quad_trial_ = quad_test_;
				reset_quadrature_ = false;
			} else {
				assert((static_cast<bool>(quad_trial_)));
				assert((static_cast<bool>(quad_test_)));
			}
		}

		class FEInitializer {
		public:

			FEInitializer(LibMeshAssemblyContext &ctx)
			: ctx(ctx)
			{}

			template<class Any>
			inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }
			

			void init_phi(const LibMeshFunctionSpace &s)
			{
				ctx.fe()[s.subspace_id()]->get_phi();
			}

			void init_phi(const ProductFunctionSpace<LibMeshFunctionSpace> &s)
			{
				s.each([&](const int, const LibMeshFunctionSpace &space) {
					ctx.fe()[space.subspace_id()]->get_phi();
				});
			}


			void init_JxW(const LibMeshFunctionSpace &s)
			{
				if(s.subspace_id() == 0) {
					//FIXME
					ctx.fe()[s.subspace_id()]->get_JxW();
				}
			}

			void init_JxW(const ProductFunctionSpace<LibMeshFunctionSpace> &s)
			{
				s.each([&](const int, const LibMeshFunctionSpace &space) {
					ctx.fe()[space.subspace_id()]->get_JxW();
				});
			}

			void init_xyz()
			{
				//FIXME
				ctx.fe()[0]->get_xyz();
			}

			template<class T>
			inline int visit(const TrialFunction<T> &expr)
			{
				init_phi(*expr.space_ptr());
				return TRAVERSE_CONTINUE;
			}

			template<class T>
			inline int visit(const TestFunction<T> &expr)
			{
				init_JxW(*expr.space_ptr());
				init_phi(*expr.space_ptr());
				return TRAVERSE_CONTINUE;
			}

			void init_dphi(const LibMeshFunctionSpace &s)
			{
				ctx.fe()[s.subspace_id()]->get_dphi();
			}

			void init_dphi(const ProductFunctionSpace<LibMeshFunctionSpace> &s)
			{
				s.each([&](const int, const LibMeshFunctionSpace &space) {
					ctx.fe()[space.subspace_id()]->get_dphi();
				});
			}

			//Gradient
			template<template<class> class Function>
			inline int visit(const Gradient<Function<ProductFunctionSpace<LibMeshFunctionSpace>>> &expr)
			{
				init_dphi(*expr.expr().space_ptr());
				return TRAVERSE_CONTINUE;
			}

			template<template<class> class Function>
			inline int visit(const Gradient<Function<LibMeshFunctionSpace>> &expr)
			{
				init_dphi(*expr.expr().space_ptr());
				return TRAVERSE_CONTINUE;
			}

			template<class C, template<class> class Function>
			inline int visit(const Gradient<Interpolate<C, Function<LibMeshFunctionSpace>> > &expr)
			{
				init_dphi(*expr.expr().fun().space_ptr());
				return TRAVERSE_SKIP_SUBTREE;
			}

			template<class C, template<class> class Function>
			inline int visit(const Gradient<Interpolate<C, Function<ProductFunctionSpace<LibMeshFunctionSpace>>> > &expr)
			{
				init_dphi(*expr.expr().fun().space_ptr());
				return TRAVERSE_SKIP_SUBTREE;
			}

			template<class C, template<class> class Function>
			inline int visit(const Interpolate<C, Function<LibMeshFunctionSpace>> &expr)
			{
				init_phi(*expr.fun().space_ptr());
				return TRAVERSE_CONTINUE;
			}

			template<class C, template<class> class Function>
			inline int visit(const Interpolate<C, Function<ProductFunctionSpace<LibMeshFunctionSpace>>> &expr)
			{
				init_phi(*expr.fun().space_ptr());
				return TRAVERSE_CONTINUE;
			}


			//Divergence
			template<template<class> class Function>
			inline int visit(const Divergence<Function<ProductFunctionSpace<LibMeshFunctionSpace>>> &expr)
			{
				init_dphi(*expr.expr().space_ptr());
				return TRAVERSE_CONTINUE;
			}

			template<template<class> class Function>
			inline int visit(const Divergence<Function<LibMeshFunctionSpace>> &expr)
			{
				init_dphi(*expr.expr().space_ptr());
				return TRAVERSE_CONTINUE;
			}

			//Curl
			template<template<class> class Function>
			inline int visit(const Curl<Function<ProductFunctionSpace<LibMeshFunctionSpace>>> &expr)
			{
				init_dphi(*expr.expr().space_ptr());
				return TRAVERSE_CONTINUE;
			}

			template<template<class> class Function>
			inline int visit(const Curl<Function<LibMeshFunctionSpace>> &expr)
			{
				init_dphi(*expr.expr().space_ptr());
				return TRAVERSE_CONTINUE;
			}

			template<class Out, class F>
			inline int visit(const ContextFunction<Out, F> &expr)
			{
				init_xyz();
				return TRAVERSE_CONTINUE;
			}

			template<class Expr>
			void apply(const Expr &expr)
			{
				traverse(expr, *this);
			}

			LibMeshAssemblyContext &ctx;
		};
	};

	template<>
	class AssemblyContext<LIBMESH_TAG> : public LibMeshAssemblyContext {};
}

#endif //UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP

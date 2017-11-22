#ifndef UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP
#define UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP 

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_Types.hpp"
#include "utopia_libmesh_FunctionSpace.hpp"
#include "utopia_libmesh_FunctionalTraits.hpp"
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

		inline std::shared_ptr<libMesh::QBase> quad_test()
		{
			return quad_test_;
		}

		inline std::shared_ptr<libMesh::QBase> quad_trial()
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
			std::cout << "init_bilinear: quadrature_order: " << quadrature_order_ << std::endl;
		}

		template<class Expr>
		void init_linear(const Expr &expr) 
		{	
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function"  );	
			init_fe_from(expr);
			std::cout << "init_linear: quadrature_order: " << quadrature_order_ << std::endl;
		}


		// static void collect_dof_indices(const libMesh::Elem * elem_ptr,
		// 							    const ProductFunctionSpace<LibMeshFunctionSpace> &space,
		// 							    std::vector<libMesh::dof_id_type> &prod_indices)
		// {
		// 	std::vector<libMesh::dof_id_type> indices;
		// 	prod_indices.clear();

		// 	const auto &sub_0 = space.subspace(0);
		// 	const auto &dof_map  = sub_0.dof_map();

		// 	space.each([&](const int sub_index, const LibMeshFunctionSpace &space) {
		// 		dof_map.dof_indices(elem_ptr, indices, space.subspace_id());
		// 		prod_indices.insert(prod_indices.end(), indices.begin(), indices.end());
		// 	});
		// }

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
			init_offsets(expr);
			
			const int dim = space_ptr->mesh().mesh_dimension();
			const libMesh::Elem * elem = space_ptr->mesh().elem(current_element_);
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

		// template<class Expr>
		// void init_test_fe(const Expr &expr)
		// {
		// 	static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function" );
		// 	auto test_space_ptr = test_space<LibMeshFunctionSpace>(expr);



		// 	if(test_space_ptr) {
		// 		test_space_ptr->initialize();
		// 		quadrature_order_ = functional_order(expr, *this);
		// 		not_null_test_ = test_space_ptr->subspace_id();
				
		// 		test_.resize(test_space_ptr->equation_system().n_vars());

		// 		const int dim = test_space_ptr->mesh().mesh_dimension();
		// 		const libMesh::Elem * elem = test_space_ptr->mesh().elem(current_element_);
		// 		set_up_quadrature(dim, quadrature_order_, elem);

		// 		auto test_fe = libMesh::FEBase::build(dim, test_space_ptr->type());
		// 		test_fe->attach_quadrature_rule(quad_test().get());
		// 		init_test_fe_flags(expr, *test_fe);
		// 		test_fe->reinit(elem);
				
		// 		test_[test_space_ptr->subspace_id()] = std::move(test_fe);
		// 		block_id_ = elem->subdomain_id();

		// 		test_space_ptr->dof_map().dof_indices(elem, test_dof_indices, test_space_ptr->subspace_id());
		// 		// n_local_test_dofs = test_space_ptr->dof_map().n_local_dofs();

		// 	} else {
		// 		auto prod_test_space_ptr = test_space<ProductFunctionSpace<LibMeshFunctionSpace>>(expr);
		// 		assert(prod_test_space_ptr.get());


		// 		prod_test_space_ptr->each([](const int, LibMeshFunctionSpace &subspace) {
		// 			subspace.initialize();
		// 		});

		// 		quadrature_order_ = functional_order(expr, *this);

		// 		//use sub0 for init geometric info
		// 		const auto &sub_0 = prod_test_space_ptr->subspace(0);
		// 		not_null_test_ = sub_0.subspace_id();
		// 		test_.resize(sub_0.equation_system().n_vars());
		// 		const int dim = sub_0.mesh().mesh_dimension();
		// 		const libMesh::Elem * elem = sub_0.mesh().elem(current_element_);
		// 		set_up_quadrature(dim, quadrature_order_, elem);
		// 		block_id_ = elem->subdomain_id();

		// 		prod_test_space_ptr->each([&](const int, LibMeshFunctionSpace &subspace) {
		// 			auto test_fe = libMesh::FEBase::build(dim, subspace.type());
		// 			test_fe->attach_quadrature_rule(quad_test().get());
		// 			init_test_fe_flags(expr, *test_fe);
		// 			test_fe->reinit(elem);
		// 			test_[subspace.subspace_id()] = std::move(test_fe);
		// 		});

		// 		collect_dof_indices(elem, *prod_test_space_ptr, test_dof_indices);
		// 		// n_local_test_dofs = sub_0.dof_map().n_local_dofs();
		// 	}


		// 	init_offsets(expr);
		// }



		// template<class Expr>
		// void init_trial_fe(const Expr &expr)
		// {
		// 	static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find trial function" );
		// 	trial_.clear();

		// 	auto trial_space_ptr = trial_space<LibMeshFunctionSpace>(expr);
		// 	auto test_space_ptr  = test_space<LibMeshFunctionSpace>(expr);

		// 	if(trial_space_ptr && trial_space_ptr == test_space_ptr) return;

		// 	if(trial_space_ptr) {
		// 		trial_space_ptr->initialize();
		// 		trial_.resize(trial_space_ptr->equation_system().n_vars());

		// 		const int dim = trial_space_ptr->mesh().mesh_dimension();
		// 		const libMesh::Elem * elem = trial_space_ptr->mesh().elem(current_element_);

		// 		auto trial_fe = libMesh::FEBase::build(dim, trial_space_ptr->type());
		// 		trial_fe->attach_quadrature_rule(quad_trial().get());
		// 		init_trial_fe_flags(expr, *trial_fe);
		// 		trial_fe->reinit(elem);
				
		// 		trial_[trial_space_ptr->subspace_id()] = std::move(trial_fe);
		// 		block_id_ = elem->subdomain_id();

		// 		trial_space_ptr->dof_map().dof_indices(elem, trial_dof_indices, trial_space_ptr->subspace_id());
		// 		// n_local_trial_dofs = trial_space_ptr->dof_map().n_local_dofs();
		// 	} else {
		// 		auto prod_trial_space_ptr = trial_space<ProductFunctionSpace<LibMeshFunctionSpace>>(expr);
		// 		auto prod_test_space_ptr  = test_space<ProductFunctionSpace<LibMeshFunctionSpace>>(expr);

		// 		if(prod_trial_space_ptr == prod_test_space_ptr) return;

		// 		assert(prod_trial_space_ptr.get());

		// 		prod_trial_space_ptr->each([](const int, LibMeshFunctionSpace &subspace) {
		// 			subspace.initialize();
		// 		});

		// 		quadrature_order_ = functional_order(expr, *this);

		// 		//use sub0 for init geometric info
		// 		const auto &sub_0 = prod_trial_space_ptr->subspace(0);
		// 		trial_.resize(sub_0.equation_system().n_vars());
		// 		const int dim = sub_0.mesh().mesh_dimension();
		// 		const libMesh::Elem * elem = sub_0.mesh().elem(current_element_);
		// 		set_up_quadrature(dim, quadrature_order_, elem);
		// 		block_id_ = elem->subdomain_id();

		// 		prod_trial_space_ptr->each([&](const int, LibMeshFunctionSpace &subspace) {
		// 			auto trial_fe = libMesh::FEBase::build(dim, subspace.type());
		// 			trial_fe->attach_quadrature_rule(quad_trial().get());
		// 			init_trial_fe_flags(expr, *trial_fe);
		// 			trial_fe->reinit(elem);
		// 			trial_[subspace.subspace_id()] = std::move(trial_fe);
		// 		});


		// 		collect_dof_indices(elem, *prod_trial_space_ptr, trial_dof_indices);
		// 		// n_local_trial_dofs = sub_0.dof_map().n_local_dofs();
		// 	} 
		// }

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

		LibMeshAssemblyContext()
		: current_element_(0), quadrature_order_(2), block_id_(0), reset_quadrature_(true)
		{}

		//x basis function
		std::vector<int> offset;
	private:
		long quadrature_order_;
		long current_element_;
		int block_id_;
		bool reset_quadrature_;


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

		// inline std::size_t n_test_functions() const
		// {
		// 	std::size_t ret = 0;
		// 	for(auto &t_ptr : test()) {
		// 		if(t_ptr) ret += t_ptr->n_shape_functions();
		// 	}

		// 	return ret;
		// }

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

		template<class Expr>
		void init_test_fe_flags(const Expr &, FE &test_fe) {
			//FIXME detect if necessary to compute
			test_fe.get_xyz();
			test_fe.get_JxW();

			if(IsSubTree<Gradient<TestFunction<utopia::Any>>, Expr>::value) {
				test_fe.get_dphi();
			}

			if(IsSubTree<TestFunction<utopia::Any>, Expr>::value) {
				test_fe.get_phi();
			}

			if(IsSubTree<Divergence<TestFunction<utopia::Any>>, Expr>::value) {
				// test_fe.get_div_phi();
				test_fe.get_dphi();
			}

			if(IsSubTree<Curl<TestFunction<utopia::Any>>, Expr>::value) {
				// test_fe.get_curl_phi();
				test_fe.get_dphi();
			}
		}

		template<class Expr>
		void init_trial_fe_flags(const Expr &, FE &trial_fe) {
			//FIXME detect if necessary to compute
			// trial_fe.get_xyz();

			if(IsSubTree<Gradient<TrialFunction<utopia::Any>>, Expr>::value) {
				trial_fe.get_dphi();
			}

			if(IsSubTree<TrialFunction<utopia::Any>, Expr>::value) {
				trial_fe.get_phi();
			}

			if(IsSubTree<Divergence<TrialFunction<utopia::Any>>, Expr>::value) {
				// trial_fe.get_div_phi();
				trial_fe.get_dphi();
			}

			if(IsSubTree<Curl<TrialFunction<utopia::Any>>, Expr>::value) {
				// trial_fe.get_curl_phi();
				trial_fe.get_dphi();
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

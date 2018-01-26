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
#include "utopia_libmesh_AssemblyValues.hpp"
#include "utopia_libmesh_TreeNavigator.hpp"

#include "utopia_fe_core.hpp"

#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/fe_interface.h"

namespace utopia {

	class HasSurfaceIntegral {
	public:

		HasSurfaceIntegral()
		: has_surface_integral(false)
		{}

		template<class Any>
		inline constexpr static int visit(const Any &) { return TRAVERSE_CONTINUE; }
		
		template<class Expr>
		inline int visit(const Integral<Expr> &expr) 
		{
			if(expr.is_surface()) {
				has_surface_integral = true;
				return TRAVERSE_STOP;
			} else {
				return TRAVERSE_CONTINUE;
			}
		}

		template<class Expr>
		void apply(const Expr &expr)
		{
			traverse(expr, *this);
		}

		bool has_surface_integral;
	};

	template<class Expr>
	bool has_surface_integral(const Expr &expr)
	{
		HasSurfaceIntegral action;
		action.apply(expr);
		return action.has_surface_integral;
	}

	class LibMeshAssemblyContext {
	public:
		typedef utopia::Traits<LibMeshFunctionSpace> TraitsT;
		typedef TraitsT::FE FE; 
		typedef TraitsT::Matrix Matrix;
		typedef TraitsT::Vector Vector;
		typedef TraitsT::DXType DXType;

		inline std::vector< std::unique_ptr<FE> > &fe()		
		{
			return active_values().fe();
		}

		inline const std::vector< std::unique_ptr<FE> > &fe() const
		{
			return active_values().fe();
		}

		inline std::vector< std::unique_ptr<FE> > &test()		
		{
			return active_values().test();
		}

		inline std::vector< std::unique_ptr<FE> > &trial()
		{
			return active_values().trial();
		}

		inline const std::vector< std::unique_ptr<FE> > &test()	const
		{
			return active_values().test();
		}

		inline const std::vector< std::unique_ptr<FE> > &trial() const
		{
			return active_values().trial();
		}

		inline std::shared_ptr<libMesh::QBase> quad_test() const
		{
			return active_values().quad_test();
		}

		inline std::shared_ptr<libMesh::QBase> quad_trial() const
		{
			return active_values().quad_trial();
		}

		template<class Expr>
		void init_bilinear(const Expr &expr) 
		{		
			static_assert( (IsSubTree<TrialFunction<utopia::Any>, Expr>::value), "could not find trial function" );
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function"  );
			init_fe_from(expr);
		}

		template<class Expr>
		void init_linear(const Expr &expr) 
		{	
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function"  );	
			init_fe_from(expr);
		}

		template<class Expr>
		void init_offsets(const Expr &expr) 
		{
			auto space_ptr = find_any_space(expr);

			assert((static_cast<bool>(space_ptr)));

			const std::size_t n_vars = space_ptr->equation_system().n_vars();
			offset.resize(n_vars + 1);
			offset[0] = 0;

			const libMesh::Elem * elem = space_ptr->mesh().elem(active_values().current_element());

			const libMesh::ElemType type = elem->type();
			const unsigned int sys_num   = space_ptr->dof_map().sys_number();
			const unsigned int dim       = elem->dim();

			for(std::size_t i = 0; i < n_vars; ++i) {
				const auto & var = space_ptr->dof_map().variable(i);
				libMesh::FEType fe_type = var.type();
				offset[i+1] = offset[i] + libMesh::FEInterface::n_shape_functions(dim, fe_type, type);
			}
		}

		static inline std::size_t n_boundary_sides(const libMesh::Elem * elem)
		{
			std::size_t ret = 0;
			for(std::size_t side = 0; side < elem->n_sides(); ++side) {
				if((elem->neighbor_ptr(side) != libmesh_nullptr)) { continue; }
				++ret;
			}

			return ret;
		}

		template<class Expr>
		void init_all_side_fe_from(const Expr &expr)
		{
			if(!has_surface_integral(expr)) {
				return;
			}

			auto space_ptr = find_any_space(expr);
			const libMesh::Elem * elem = space_ptr->mesh().elem(active_values().current_element());

			auto n = n_boundary_sides(elem);

			surface_values_.resize(n);

			std::size_t index = 0;
			for(std::size_t side = 0; side < elem->n_sides(); ++side) {
				if((elem->neighbor_ptr(side) != libmesh_nullptr)) { continue; }

				auto &values_ptr = surface_values_[index++];
				
				if(!values_ptr) {
					values_ptr = std::make_shared<LibMeshAssemblyValues>();
				}

				values_ptr->set_current_element(active_values().current_element());
				values_ptr->init_side_fe_from(expr, side);
			}

		}

		template<class Expr>
		void init_fe_from(const Expr &expr)
		{
			active_values().init_fe_from(expr);
			init_offsets(expr);

			init_all_side_fe_from(expr);
		}

		void init_tensor(Vector &v, const bool reset);
		void init_tensor(Matrix &v, const bool reset);

		const DXType &dx() const
		{
			return test()[0]->get_JxW();
		}

		inline int block_id() const
		{
			return active_values().block_id();
		}

		inline void set_current_element(const long current_element)
		{
			active_values().set_current_element(current_element);
		}

		inline long current_element() const
		{
			return active_values().current_element();
		}

		void surface_integral_begin() {}

		void surface_integral_end()
		{
			active_values_ = volume_values_;
		}

		void set_side(const std::size_t i)
		{
			assert(surface_values_[i]);
			active_values_ = surface_values_[i];
		}

		std::size_t n_sides()
		{
			return surface_values_.size();
		}

		LibMeshAssemblyContext()
		: has_assembled_(false)
		{
			active_values_ = std::make_shared<LibMeshAssemblyValues>();
		}

		//x basis function
		std::vector<int> offset;

		inline void set_has_assembled(const bool val)
		{
			has_assembled_ = val;
		}

		inline bool has_assembled() const
		{
			return has_assembled_;
		}

	private:
		std::shared_ptr<libMesh::QBase> quad_trial_;
		std::shared_ptr<libMesh::QBase> quad_test_;
		std::vector< std::unique_ptr<FE> > fe_;
		std::shared_ptr<LibMeshAssemblyValues> active_values_;
		std::shared_ptr<LibMeshAssemblyValues> volume_values_;
		std::vector< std::shared_ptr<LibMeshAssemblyValues> > surface_values_;
		bool has_assembled_;

		inline LibMeshAssemblyValues &active_values()
		{
			assert(active_values_);
			return *active_values_;
		}

		inline const LibMeshAssemblyValues &active_values() const
		{
			assert(active_values_);
			return *active_values_;
		}

		inline std::size_t n_shape_functions() const
		{
			return active_values().n_shape_functions();
		}
	};

	template<>
	class AssemblyContext<LIBMESH_TAG> : public LibMeshAssemblyContext {};
}

#endif //UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP

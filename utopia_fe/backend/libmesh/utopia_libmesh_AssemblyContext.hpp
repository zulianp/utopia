#ifndef UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP
#define UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP 

#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_fe_core.hpp"

namespace utopia {

	class LibMeshAssemblyContext {
	public:
		typedef utopia::LibMeshTraits::FE FE; 
		typedef utopia::LibMeshTraits::Matrix Matrix;
		typedef utopia::LibMeshTraits::Vector Vector;
		typedef utopia::LibMeshTraits::DxType DXType;

		std::vector< std::unique_ptr<FE> > trial_;
		std::vector< std::unique_ptr<FE> > test_;

		long quadrature_order_;
		long current_element_;
		int block_id_;

		inline std::vector< std::unique_ptr<FE> > &test()		
		{
			return test_;
		}

		inline std::vector< std::unique_ptr<FE> > &trial()
		{
			if(trial_.empty()) {
				//if space is symmetric
				return test();
			} else {
				return trial_;
			}
		}

		inline std::shared_ptr<libMesh::QGauss> quad_test()
		{
			return quad_test_;
		}

		inline std::shared_ptr<libMesh::QGauss> quad_trial()
		{
			return quad_trial_;
		}

		template<class Expr>
		void init_bilinear(const Expr &expr) 
		{		
			static_assert( (IsSubTree<TrialFunction<utopia::Any>, Expr>::value), "could not find trial function" );
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function" );
			//TODO
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
				test_fe.get_div_phi();
			}

			if(IsSubTree<Curl<TestFunction<utopia::Any>>, Expr>::value) {
				test_fe.get_curl_phi();
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
				trial_fe.get_div_phi();
			}

			if(IsSubTree<Curl<TrialFunction<utopia::Any>>, Expr>::value) {
				trial_fe.get_curl_phi();
			}
		}

		template<class Expr>
		void init_linear(const Expr &expr) 
		{		
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function" );

			quadrature_order_ = functional_order(expr, *this);

			auto test_space_ptr = test_space<LibMeshFunctionSpace>(expr);

			if(test_space_ptr) {
				test_.resize(test_space_ptr->equation_system().n_vars());

				const int dim = test_space_ptr->mesh().mesh_dimension();
				const libMesh::Elem * elem = test_space_ptr->mesh().elem(current_element_);
				set_up_quadrature(dim, quadrature_order_, elem);

				auto test_fe = libMesh::FEBase::build(dim, test_space_ptr->type());
				test_fe->attach_quadrature_rule(quad_test().get());
				init_test_fe_flags(expr, *test_fe);
				test_fe->reinit(elem);
				
				test_[test_space_ptr->subspace_id()] = std::move(test_fe);
				block_id_ = elem->subdomain_id();
			} else {
				auto prod_test_space_ptr = test_space<ProductFunctionSpace<LibMeshFunctionSpace>>(expr);
				assert(prod_test_space_ptr);
			}
		}

		void init_tensor(Vector &v, const bool reset);
		void init_tensor(Matrix &v, const bool reset);

		const DXType &dx() const
		{
			return test[0]->get_JxW();
		}

		AssemblyContext()
		: current_element_(0), quadrature_order_(2), block_id_(0), reset_quadrature_(true)
		{}

		inline int block_id() const
		{
			return block_id_;
		}

		bool reset_quadrature_;
		std::shared_ptr<libMesh::QGauss> quad_trial_;
		std::shared_ptr<libMesh::QGauss> quad_test_;

	private:
		void set_up_quadrature(const int dim, const int quadrature_order, libMesh::Elem * elem)
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
	};

	template<>
	class AssemblyContext<LIBMESH_TAG> : public LibMeshAssemblyContext {};
}

#endif //UTOPIA_LIBMESH_ASSEMBLY_CONTEXT_HPP

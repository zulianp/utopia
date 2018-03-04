#ifndef UTOPIA_HOMEMADE_ASSEMBLY_CONTEXT_HPP
#define UTOPIA_HOMEMADE_ASSEMBLY_CONTEXT_HPP 


#include "utopia_Base.hpp"
#include "utopia_Traits.hpp"
#include "utopia_AssemblyContext.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_Traverse.hpp"
#include "utopia_ProductFunctionSpace.hpp"


#include "utopia_homemade_FEForwardDeclarations.hpp"
#include "utopia_homemade_FE.hpp"


#include "utopia_fe_core.hpp"

namespace utopia {

	template<>
	class AssemblyContext<HOMEMADE> {
	public:

		std::vector< std::shared_ptr<FE> > trial;
		std::vector< std::shared_ptr<FE> > test;

		long quadrature_order;
		long current_element;
		int block_id_;

		template<class Expr>
		void init_bilinear(const Expr &expr) 
		{		

			static_assert( (IsSubTree<TrialFunction<utopia::Any>, Expr>::value), "could not find trial function" );
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), "could not find test function" );

			//clean up previous context
			trial.clear();
			test.clear();

			std::shared_ptr<FE> trial_fe;
			std::shared_ptr<FE> test_fe;

			auto trial_space_ptr = trial_space<HMFESpace>(expr);
			auto test_space_ptr  = test_space<HMFESpace>(expr);


			quadrature_order = functional_order(expr, *this);

			assert(functional_type(expr, *this) == POLYNOMIAL_FUNCTION);

			std::cout << "quadrature_order: " << quadrature_order << std::endl;

			if(trial_space_ptr) {
				trial_fe = std::make_shared<FE>();
				trial_fe->init(current_element, trial_space_ptr->mesh(), quadrature_order);
				trial.push_back(trial_fe);
			} else {
				//product space init
				auto prod_trial_space_ptr = trial_space<ProductFunctionSpace<HMFESpace>>(expr);
				assert(prod_trial_space_ptr);

				trial_fe = std::make_shared<FE>();
				trial_fe->init(current_element, prod_trial_space_ptr->subspace(0).mesh(), quadrature_order);

				//just copy it three times for now
				for(std::size_t i = 0; i < prod_trial_space_ptr->n_subspaces(); ++i) {
					trial.push_back(trial_fe);
				}
			}

			if(test_space_ptr) {
				if(trial_space_ptr != test_space_ptr) {
					test_fe = std::make_shared<FE>();
					test_fe->init(current_element, test_space_ptr->mesh(), quadrature_order);
					test.push_back(test_fe);
				} else if(trial_space_ptr) {
					test_fe = trial_fe;
				}
			} else {
				//product space init
				auto prod_test_space_ptr = test_space<ProductFunctionSpace<HMFESpace>>(expr);
				assert(prod_test_space_ptr);

				test_fe = std::make_shared<FE>();
				test_fe->init(current_element, prod_test_space_ptr->subspace(0).mesh(), quadrature_order);

				//just copy it three times for now
				for(std::size_t i = 0; i < prod_test_space_ptr->n_subspaces(); ++i) {
					test.push_back(test_fe);
				}
			}
		}

		template<class Expr>
		void init_linear(const Expr &expr) 
		{		
			static_assert( (IsSubTree<TestFunction<utopia::Any>,  Expr>::value), 	"could not find test function" );

			test.clear();
			
			auto test_space_ptr  = test_space<HMFESpace>(expr);
			std::shared_ptr<FE> test_fe;

			if(test_space_ptr) {
				test_fe = std::make_shared<FE>();
				test_fe->init(current_element, test_space_ptr->mesh(), quadrature_order);
				test.push_back(test_fe);
			
			} else {
				//product space init
				auto prod_test_space_ptr = test_space<ProductFunctionSpace<HMFESpace>>(expr);
				assert(prod_test_space_ptr);

				test_fe = std::make_shared<FE>();
				test_fe->init(current_element, prod_test_space_ptr->subspace(0).mesh(), quadrature_order);

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
		: current_element(0), quadrature_order(2), block_id_(0)
		{}

		inline int block_id() const
		{
			return block_id_;
		}
	};
}

#endif //UTOPIA_HOMEMADE_ASSEMBLY_CONTEXT_HPP

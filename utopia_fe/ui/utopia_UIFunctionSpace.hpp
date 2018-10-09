#ifndef UTOPIA_UI_FUNCTION_SPACE_HPP
#define UTOPIA_UI_FUNCTION_SPACE_HPP

#include "utopia_ui.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_UIMesh.hpp"

#include "libmesh/mesh_refinement.h"

#include <memory>

namespace utopia {

	template<class FunctionSpace>
	class UIFunctionSpace final : public Configurable {
	public:
		UIFunctionSpace()
		{}

		void read(Input &is) override {}

	private:
		std::shared_ptr<FunctionSpace> space_;
	};


	template<>
	class UIFunctionSpace<LibMeshFunctionSpace> final : public Configurable {
	public:
		UIFunctionSpace(const std::shared_ptr<UIMesh<libMesh::DistributedMesh>> &mesh)
		: mesh_(mesh)
		{}

		void read(Input &is) override {
			int n_vars = 0;

			std::vector<std::string> var_names;
			std::vector<int> var_orders;

			std::string system_name = "main";

			is.read("variables", [&](Input &sub_is_mid) {
				sub_is_mid.read_all([&](Input &sub_is) {
					std::string var_name = "u_" + std::to_string(n_vars);
					int var_order = 1;

					sub_is.read("name", var_name);
					sub_is.read("order", var_order);
					
					var_names.push_back(var_name);
					var_orders.push_back(var_order);

					++n_vars;
				});
			});

			if(n_vars == 0) {
				var_names.push_back("u");
				var_orders.push_back(1);
				n_vars = 1;
			}

			is.read("system-name", system_name);

			auto equation_systems = std::make_shared<libMesh::EquationSystems>(mesh_->mesh());
			auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>(system_name);

			space_ = std::make_shared<ProductFunctionSpace<LibMeshFunctionSpace>>();

			for(int i = 0; i < n_vars; ++i) {
				auto ss = std::make_shared<LibMeshFunctionSpace>(equation_systems, libMesh::LAGRANGE, libMesh::Order(var_orders[i]), var_names[i]);
				space_->add_subspace(ss);
			}

			is.read("boundary-conditions", [this](Input &is) {
			    is.read_all([this](Input &is) {
			        int side_set = 0;
			        
			        is.read("side", side_set);
			        
			        double value = 0;
			        is.read("value", value);

			        int var_num = 0;

			        is.read("var", var_num);

			        auto u = trial(space_->subspace(var_num));

			        init_constraints(
			        	constraints(
			        		boundary_conditions(u == coeff(value), {side_set})
			        		)
			        	);

			    });
			});

			for(int i = 0; i < n_vars; ++i) {
				space_->subspace(i).initialize();
			}
		}

		inline ProductFunctionSpace<LibMeshFunctionSpace> &space()
		{
			assert(space_);
			return *space_;
		}

		inline LibMeshFunctionSpace &subspace(int index)
		{
			return space_->subspace(index);
		}

		inline std::shared_ptr<LibMeshFunctionSpace> subspace_ptr(int index)
		{
			return space_->subspace_ptr(index);
		}

	private:
		std::shared_ptr<UIMesh<libMesh::DistributedMesh>> mesh_;
		std::shared_ptr<ProductFunctionSpace<LibMeshFunctionSpace>> space_;
	};
}


#endif //UTOPIA_UI_FUNCTION_SPACE_HPP

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
		UIFunctionSpace(
			const std::shared_ptr<UIMesh<libMesh::DistributedMesh>> &mesh,
			const std::shared_ptr<libMesh::EquationSystems> equation_systems = nullptr)
		: mesh_(mesh), equation_systems_(equation_systems)
		{}

		void read(Input &is) override {
			int n_vars = 0;

			std::vector<std::string> var_names;
			std::vector<int> var_orders;
			std::vector<std::string> fe_families;

			std::string system_name = "main";
			std::string fe_family = "LAGRANGE";

			is.get("variables", [&](Input &sub_is_mid) {
				sub_is_mid.get_all([&](Input &sub_is) {
					std::string var_name = "u_" + std::to_string(n_vars);
					int var_order = 1;
					std::string var_fe_family = fe_family;

					sub_is.get("name", var_name);
					sub_is.get("order", var_order);
					sub_is.get("fe-family", var_fe_family);
					
					var_names.push_back(var_name);
					var_orders.push_back(var_order);
					fe_families.push_back(var_fe_family);

					++n_vars;
				});
			});

			if(n_vars == 0) {
				var_names.push_back("u");
				var_orders.push_back(1);
				fe_families.push_back(fe_family);
				n_vars = 1;
			}

			is.get("system-name", system_name);

			if(!equation_systems_) {
				equation_systems_ = std::make_shared<libMesh::EquationSystems>(mesh_->mesh());
			}

			auto &sys = equation_systems_->add_system<libMesh::LinearImplicitSystem>(system_name);
			space_ = std::make_shared<ProductFunctionSpace<LibMeshFunctionSpace>>();

			for(int i = 0; i < n_vars; ++i) {
				auto ss = std::make_shared<LibMeshFunctionSpace>(
					equation_systems_,
					libMesh::Utility::string_to_enum<libMesh::FEFamily>(fe_families[i]),
					libMesh::Order(var_orders[i]),
					var_names[i],
					sys.number()
				);
				
				space_->add_subspace(ss);
			}

			is.get("boundary-conditions", [this](Input &is) {
			    is.get_all([this](Input &is) {
			        int side_set = 0;
			        
			        is.get("side", side_set);
			        
			        double value = 0;
			        is.get("value", value);

			        int var_num = 0;

			        is.get("var", var_num);

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

		inline void set_space(const std::shared_ptr<ProductFunctionSpace<LibMeshFunctionSpace> > &space)
		{
			space_ = space;
		}

		inline LibMeshFunctionSpace &subspace(int index)
		{
			return space_->subspace(index);
		}

		inline std::shared_ptr<LibMeshFunctionSpace> subspace_ptr(int index)
		{
			return space_->subspace_ptr(index);
		}

		std::shared_ptr<UIMesh<libMesh::DistributedMesh>> mesh() 
		{
			return mesh_;
		}

		inline bool initialized() const
		{
			return static_cast<bool>(space_);
		}

	private:
		std::shared_ptr<UIMesh<libMesh::DistributedMesh>> mesh_;
		std::shared_ptr<libMesh::EquationSystems> equation_systems_;
		std::shared_ptr<ProductFunctionSpace<LibMeshFunctionSpace>> space_;
	};
}


#endif //UTOPIA_UI_FUNCTION_SPACE_HPP

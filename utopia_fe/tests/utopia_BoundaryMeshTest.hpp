#ifndef UTOPIA_BOUNDARY_MESH_TEST_HPP
#define UTOPIA_BOUNDARY_MESH_TEST_HPP

#include "utopia_fe_base.hpp"
#include "utopia_FETest.hpp"
#include <string>


namespace utopia {

	class BoundaryMeshTest final : public FETest {
	public:
		void run(Input &in) override;

		inline static std::string command()
		{
			return "bit";
		}
	};

}


#endif //UTOPIA_BOUNDARY_MESH_TEST_HPP

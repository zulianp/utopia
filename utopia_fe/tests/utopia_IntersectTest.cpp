#include "utopia_IntersectTest.hpp"
#include "utopia_libmesh.hpp"
#include "utopia_Intersect.hpp"


namespace utopia {
	void run_intersect_test(libMesh::LibMeshInit &init)
	{
		Polygon3 poly =
		{
			{ 
				{ 0.4, 0.5, 0.3 },
				{ 0.4, 0.5, 0.2 },
				{ 0.3, 0.5, 0.2 },
				{ 0.3, 0.5, 0.3 }
			}
		}; 

		HPolyhedron3 h =
		{
			{
				{ //0
					{ 0, -1, 0},
					-0.33333333333333331
				},
				{ //1
					{ 1, 0, 0},
					0.33333333333333331
				},
				{ //2
					{ 0, 1, 0},
					0.66666666666666663
				},
				{ //3
					{ -1, 0, 0},
					0
				},
				{ //4
					{ 0, -0, -1},
					0
				},
				{ //5
					{ 0, 0, 1},
					0.33333333333333331
				}
			}
		};

		// Polygon3 poly =
		// {
		// 	{ 
		// 		{ 0.0, 0.0, 0.0 },
		// 		{ 1.0, 0.0, 0.0 },
		// 		{ 1.0, 1.0, 0.0 },
		// 		{ 0.0, 1.0, 0.0 }
		// 	}
		// }; 

		// auto s = 1./std::sqrt(3);

		// HPolyhedron3 h =
		// {
		// 	{
		// 		{ //0
		// 			{ 1., 0., 0.},
		// 			0.8
		// 		}
		// 		,
		// 		{ //1
		// 			{ 0, 1., 0.},
		// 			0.8
		// 		},
		// 		{
		// 			{ s, s, s },
		// 			0.5
		// 		}


		// 	}
		// };


		Polygon3 isect;

		// poly.plot("test");

		bool ok = intersect(poly, h, isect, 1e-10);
		assert(ok);


		Polygon3 oracle =
		{
			{
				{ 0.33333333333333331, 0.5, 0.20000000000000001 },
				{ 0.29999999999999999, 0.5, 0.20000000000000001 },
				{ 0.29999999999999999, 0.5, 0.29999999999999999 },
				{ 0.33333333333333331, 0.5, 0.29999999999999999 }
			}
		};

		utopia_test_assert(oracle.equals(isect, 1e-10));
		
		// isect.plot("isect");

	}
}
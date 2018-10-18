#include "utopia_libmesh.hpp"
#include "utopia_IntersectTest.hpp"
#include "utopia_Intersect.hpp"


namespace utopia {
	void intersect_hex_with_polygon_test()
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

		Polygon3 isect;

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
	}

	void intersect_tet_with_polygon_test_1()
	{
		Polygon3 poly =
		{
			{ 
				{0.5, 0.3, 0.},
				{0.5, 0.3, 0.1},
				{0.5, 0.2, 0.1},
				{0.5, 0.2, 0.}
			}
		}; 

		HPolyhedron3 h =
		{
			{
				{
					{0.70710678118654746,  0.70710678118654746,  0},
					0.70710678118654746
				}
				,
				{
					{-0.70710678118654746,  0.70710678118654746,  0},
					0.
				}
				,
				{
					{0,  -0.70710678118654746,  0.70710678118654746},
					0.
				}
				,
				{
					{0, 0, 1.},
					0.
				}
			}
		};

		Polygon3 isect;

		bool ok = intersect(poly, h, isect, 1e-10);
		assert(!ok);

		//no intersection

		// poly.plot("poly");
		// isect.plot("isect");
	}

	void intersect_tet_with_polygon_test_2()
	{
		Polygon3 poly =
		{
			{ 
				{0.5, 0.3, 0.},
				{0.5, 0.3, 0.1},
				{0.5, 0.2, 0.1},
				{0.5, 0.2, 0.}
			}
		}; 

		HPolyhedron3 h =
		{
			{
				{
					{0.70710678118654746,  0.70710678118654746,  0},
					0.70710678118654746
				}
				,
				{
					{-0.70710678118654746,  0.70710678118654746,  0},
					0.
				}
				,
				{
					{0,  -0.70710678118654746,  0.70710678118654746},
					0.
				}
				,
				{
					{0, 0, 1.},
					0.05
				}
			}
		};

		Polygon3 isect;

		bool ok = intersect(poly, h, isect, 1e-10);
		assert(ok);

		//add oracle

		// poly.plot("poly");
		// isect.plot("isect");
	}


	void intersect_tet_with_polygon_test_3()
	{
		Polygon3 poly =
		{
			{ 
				// { 0.5, 0.099999999999999978, 0.099999999999999978 },
				// { 0.5, 0.099999999999999978, 0.19999999999999996 },
				// { 0.5, 0, 0.19999999999999996 },
				// { 0.5, 0, 0.099999999999999978 }
				{ 0.5, 0.1, 0.1 },
				{ 0.5, 0.1, 0.2 },
				{ 0.5, 0.,  0.2 },
				{ 0.5, 0.,  0.1 }
			}
		}; 

		HPolyhedron3 h =
		{
			{
				{
					{0.70710678118654746, 0.70710678118654746, 0.},
				    0.70710678118654746
				 }
				,
				 {
					{-0.70710678118654746, 0.70710678118654746, 0.},
				    0.
				 }
				,
				 {
					{0., -0.70710678118654746, 0.70710678118654746},
				    0.
				 }
				,
				 {
					{-0., -0., -1.},
				    -0.
				 }
			}
		};

		Polygon3 isect;

		bool ok = intersect(poly, h, isect, 1e-10);
		assert(!ok);
		//add oracle

		poly.plot("poly");
		isect.plot("isect");

		assert(ok);
	}

	void run_intersect_test(libMesh::LibMeshInit &init)
	{
		intersect_hex_with_polygon_test();
		intersect_tet_with_polygon_test_1();
		intersect_tet_with_polygon_test_2();
		intersect_tet_with_polygon_test_3();
		
		// isect.plot("isect");

	}
}

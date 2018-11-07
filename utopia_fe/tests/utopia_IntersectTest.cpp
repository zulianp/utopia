#include "utopia_libmesh.hpp"
#include "utopia_IntersectTest.hpp"
#include "utopia_Intersect.hpp"
#include "utopia_LibMeshShape.hpp"

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

		// assert(ok);
	}

	void intersect_convex_polyhedra_test()
	{
		Polyhedron poly1, poly2;


		poly1.n_elements = 6;
		poly1.n_nodes = 8;
		poly1.n_dims = 3;

		std::vector<int> el_ptr =
		{
		    0, 	//[0]
		    4, 	//[1]
		    8, 	//[2]
		    12, //[3]
		    16, //[4]
		    20, //[5]
		    24 	//[6]
		};

		std::vector<int> el_index =
		{
		    0, 	//[0]
		    1, 	//[1]
		    5, 	//[2]
		    4, 	//[3]
		    1, 	//[4]
		    2, 	//[5]
		    6, 	//[6]
		    5, 	//[7]
		    3, 	//[8]
		    7, 	//[9]
		    6, 	//[10]
		    2, 	//[11]
		    0, 	//[12]
		    4, 	//[13]
		    7, 	//[14]
		    3, 	//[15]
		    2, 	//[16]
		    1, 	//[17]
		    0, 	//[18]
		    3, 	//[19]
		    6, 	//[20]
		    7, 	//[21]
		    4, 	//[22]
		    5 	//[23]
		};

		std::vector<double>  points =
		{
		    -0.30000000000000004, 	//[0]
		    -0.30000000000000004, 	//[1]
		    0.30000000000000004, 	//[2]
		    -0.19999999999999996, 	//[3]
		    -0.30000000000000004, 	//[4]
		    0.30000000000000004, 	//[5]
		    -0.19999999999999996, 	//[6]
		    -0.19999999999999996, 	//[7]
		    0.30000000000000004, 	//[8]
		    -0.30000000000000004, 	//[9]
		    -0.19999999999999996, 	//[10]
		    0.30000000000000004, 	//[11]
		    -0.30000000000000004, 	//[12]
		    -0.30000000000000004, 	//[13]
		    0.39999999999999991, 	//[14]
		    -0.19999999999999996, 	//[15]
		    -0.30000000000000004, 	//[16]
		    0.39999999999999991, 	//[17]
		    -0.19999999999999996, 	//[18]
		    -0.19999999999999996, 	//[19]
		    0.39999999999999991, 	//[20]
		    -0.30000000000000004, 	//[21]
		    -0.19999999999999996, 	//[22]
		    0.39999999999999991 	//[23]
		};

		poly1.type = 3;

		std::copy(std::begin(el_ptr),   std::end(el_ptr),   std::begin(poly1.el_ptr));
		std::copy(std::begin(el_index), std::end(el_index), std::begin(poly1.el_index));
		std::copy(std::begin(points),   std::end(points),   std::begin(poly1.points));

		/////////////////////////////////////////

		poly2.n_elements = 6;
		poly2.n_nodes = 8;
		poly2.n_dims = 3;

		el_ptr = {
		  0,		//[0]
		  4,		//[1]
		  8,		//[2]
		  12,		//[3]
		  16,		//[4]
		  20,		//[5]
		  24		//[6]
		};

		el_index = {
		  0,		//[0]
		  1,		//[1]
		  5,		//[2]
		  4,		//[3]
		  1,		//[4]
		  2,		//[5]
		  6,		//[6]
		  5,		//[7]
		  3,		//[8]
		  7,		//[9]
		  6,		//[10]
		  2,		//[11]
		  0,		//[12]
		  4,		//[13]
		  7,		//[14]
		  3,		//[15]
		  2,		//[16]
		  1,		//[17]
		  0,		//[18]
		  3,		//[19]
		  6,		//[20]
		  7,		//[21]
		  4,		//[22]
		  5		    //[23]
		};

		points = {
		  -0.22651513583561636,		//[0]
		  -0.22651513583561636,		//[1]
		  0.29887968534141401,		//[2]
		  -0.15241447521870757,		//[3]
		  -0.22651513583561642,		//[4]
		  0.300986262333859,		//[5]
		  -0.15241447521870757,		//[6]
		  -0.15241447521870757,		//[7]
		  0.30309283932630399,		//[8]
		  -0.22651513583561636,		//[9]
		  -0.15241447521870757,		//[10]
		  0.300986262333859,		//[11]
		  -0.26504358749434054,		//[12]
		  -0.26504358749434054,		//[13]
		  0.35174008807660145,		//[14]
		  -0.17834708691207959,		//[15]
		  -0.26504358749434054,		//[16]
		  0.35421713095038032,		//[17]
		  -0.17834708691207959,		//[18]
		  -0.17834708691207959,		//[19]
		  0.35669417382415919,		//[20]
		  -0.26504358749434054,		//[21]
		  -0.17834708691207959,		//[22]
		  0.35421713095038032		//[23]
		};

		poly2.type = 3;

		std::copy(std::begin(el_ptr),   std::end(el_ptr),   std::begin(poly2.el_ptr));
		std::copy(std::begin(el_index), std::end(el_index), std::begin(poly2.el_index));
		std::copy(std::begin(points),   std::end(points),   std::begin(poly2.points));


		// std::cout << "-----------------------------\n";
		// Intersector::p_mesh_print(&poly1);
		// std::cout << "-----------------------------\n";
		// Intersector::p_mesh_print(&poly2);
		// std::cout << "-----------------------------\n";

		Polyhedron isect;
		bool ok = Intersector::intersect_convex_polyhedra(poly1, poly2, &isect);
	}

	void intesect_ray_elem_tet(libMesh::Parallel::Communicator &comm)
	{
		libMesh::Mesh mesh(comm);
		libMesh::MeshTools::Generation::build_square(mesh,
			1, 1,
			0., 1.,
			0., 2.,
			libMesh::QUAD8
		);

		LibMeshFunctionSpace V(mesh, libMesh::LAGRANGE, libMesh::SECOND);
		V.initialize();

		Ray<double, 2> ray = {
			{-0.5, 0.6 },
			{ 1., 0. }
		};

		// int q = 0;

		for(auto e_it = elements_begin(mesh); e_it != elements_end(mesh); ++e_it) {
			const auto &e = **e_it;

			auto n_sides = e.n_sides();

			for(std::size_t i = 0; i < n_sides; ++i) {
				if(i != 1 && i != 3) continue;

				auto side_ptr = e.build_side_ptr(i);
				LibMeshShape<double, 2> shape(*side_ptr, V.dof_map().variable_type(0));

				// std::cout << "elem: " << i << std::endl;
				// std::cout << "(" << side_ptr->node_ref(0)(0) << ", " << side_ptr->node_ref(0)(1) << ") ";
				// std::cout << "(" << side_ptr->node_ref(1)(0) << ", " << side_ptr->node_ref(1)(1) << ") " << std::endl;;

				double t = 0.;
				shape.intersect(ray, t);

				if(i == 1) {
					utopia_test_assert(approxeq(t, 1.5, 1e-8));
				} else if(i == 3) {
					utopia_test_assert(approxeq(t, 0.5, 1e-8));
				}

				// std::cout << "intersection: " << t << std::endl;
			}
		}
	}

	void run_intersect_test(libMesh::LibMeshInit &init)
	{
		UTOPIA_UNIT_TEST_BEGIN("IntersectTest");

		intesect_ray_elem_tet(init.comm());

		UTOPIA_RUN_TEST(intersect_hex_with_polygon_test);
		UTOPIA_RUN_TEST(intersect_tet_with_polygon_test_1);
		UTOPIA_RUN_TEST(intersect_tet_with_polygon_test_2);
		UTOPIA_RUN_TEST(intersect_tet_with_polygon_test_3);
		UTOPIA_RUN_TEST(intersect_convex_polyhedra_test);

		UTOPIA_UNIT_TEST_END("IntersectTest");
	}
}

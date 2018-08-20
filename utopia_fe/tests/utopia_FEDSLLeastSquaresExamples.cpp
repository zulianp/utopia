

// #include <iostream>

// #include "utopia_FEDSLLeastSquaresExamples.hpp"
// #include "utopia.hpp"

// //fe extension
// #include "utopia_fe_core.hpp"
// #include "MortarAssembler.hpp"
// #include "utopia_Socket.hpp"
// #include "utopia_ContactSimParams.hpp"
// #include "utopia_LibMeshBackend.hpp"


// using namespace utopia;
// using namespace std;
// using namespace libMesh;


// namespace utopia {

// 	LMDenseMatrix make_stress_strain_rel_tensor(const int dim, const double mu, const double lambda)
// 	{
// 		const int dim2 = dim * dim;
// 		LMDenseMatrix A = zeros(dim2, dim2);	

// 		switch(dim) {
// 			case 2:
// 			{	
// 				A.set(0, 0, 1);
// 				A.set(0, dim2-1, 1);
// 				A.set(dim2-1, 0, 1);
// 				A.set(dim2-1, dim2-1, 1);

// 				A *= -lambda / (dim * lambda + 2 * mu);
// 				A += identity(size(A));
// 				A *= 1./(2 * mu);			
// 				break;
// 			}

// 			case 3:
// 			{

// 				A.set(0, 0, 1);
// 				A.set(0, 4, 1);
// 				A.set(0, dim2-1, 1);

// 				A.set(4, 0, 1);
// 				A.set(4, 4, 1);
// 				A.set(4, dim2-1, 1);

// 				A.set(dim2-1, 0, 1);
// 				A.set(dim2-1, 4, 1);
// 				A.set(dim2-1, dim2-1, 1);

// 				A *= -lambda / (dim * lambda + 2 * mu);
// 				A += identity(size(A));
// 				A *= 1./(2 * mu);
// 				break;
// 			}

// 			default:
// 			{
// 				assert(false);
// 			}
// 		}

// 		return A;
// 	}

// 	// Example of local system on element 0
// 	void assemble_ls_local(LibMeshFEFunction &v, LibMeshVecFEFunction &u)
// 	{
// 		using namespace libMesh;

// 		u.set_element(0);
// 		v.set_element(0);

// 		//bilinear forms
// 		double c = -1.0;
// 		auto eq_11 = integral((c*c) * dot(v, v) + dot(grad(v), grad(v)));
// 		auto eq_12 = integral(c * dot(div(u), v) + dot(u, grad(v)));
// 		auto eq_21 = integral(c * dot(v, div(u)) + dot(grad(v), u));
// 		auto eq_22 = integral(dot(u, u) + dot(div(u), div(u)) + dot(curl(u), curl(u)));

// 		//linear forms
// 		auto f = coeff(1.0);
// 		auto rhs_1 = integral(c * dot(f, v));
// 		auto rhs_2 = integral(dot(f, div(u)));


// 		//assemble local block system
// 		DenseMatrix<Real> mat_11, mat_12, mat_21, mat_22;
// 		DenseVector<Real> vec_1, vec_2;

// 		block_local_assemble(
// 			//Equations
// 			eq_11, eq_12, rhs_1, 
// 			eq_21, eq_22, rhs_2, 
// 			//Matrices and vectors
// 			mat_11, mat_12, vec_1,
// 			mat_21, mat_22, vec_2
// 			);

// 		logger() << "block-system dimensions:\n";
// 		logger() << mat_11.m() << " x " << mat_11.n() << endl;
// 		logger() << mat_12.m() << " x " << mat_12.n() << endl;
// 		logger() << mat_21.m() << " x " << mat_21.n() << endl;
// 		logger() << mat_22.m() << " x " << mat_22.n() << endl;
// 	}

// 	void leastsquares_helmoholtz(LibMeshInit &init)
// 	{
// 		auto mesh = make_shared<libMesh::Mesh>(init.comm());		
// 		MeshTools::Generation::build_square (*mesh,
// 			40, 40,
// 			-1., 1.,
// 			-1., 1.,
// 			QUAD9);

// 		LibMeshFEContext<LinearImplicitSystem> context(mesh);
// 		auto Vh = fe_space(LAGRANGE, SECOND, context);
// 		auto Uh = vector_fe_space(LAGRANGE_VEC, FIRST, context);
		
// 		auto v  = fe_function(Vh);
// 		auto u  = fe_function(Uh);	
		
// 		strong_enforce( boundary_conditions(v == coeff(0.0), {0, 1, 2, 3}) );

// 		//system initialized
// 		context.equation_systems.init();

// 		const int dim = mesh->mesh_dimension();

// 		//FIXME to be moved in the assembly loop
// 		v.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		u.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));

// 		//bilinear forms
// 		double c = -100.0;
// 		double beta = 0.99;
// 		auto eq_11 = integral((c*c) * dot(v, v) + dot(grad(v), grad(v)));
// 		auto eq_12 = integral(c * dot(div(u), v) + dot(u, grad(v)));
// 		auto eq_21 = integral(c * dot(v, div(u)) + dot(grad(v), u));
// 		auto eq_22 = integral(dot(u, u) + dot(div(u), div(u)) + beta * dot(curl(u), curl(u)));

// 		//linear forms
// 		auto f = coeff(1);
// 		auto rhs_1 = integral(c * dot(f, v));
// 		auto rhs_2 = integral(dot(f, div(u)));

// 		auto &matrix = *context.system.matrix;
// 		auto &rhs 	 = *context.system.rhs;

// 		auto ass = make_assembly([&]() -> void {
// 			double t = MPI_Wtime();

// 			assemble(v, v, eq_11, matrix);
// 			assemble(u, v, eq_12, matrix);
// 			assemble(v, u, eq_21, matrix);
// 			assemble(u, u, eq_22, matrix);
// 			assemble(v,    rhs_1, rhs);
// 			assemble(u,    rhs_2, rhs);

// 			t = MPI_Wtime() - t;

// 			printf("--------------------------------\n");
// 			printf("Assembly: %g seconds\n", t);
// 			printf("--------------------------------\n");
// 		});

// 		context.system.attach_assemble_object(ass);
// 		context.equation_systems.print_info();
// 		context.equation_systems.solve();

// 		ExodusII_IO(*mesh).write_equation_systems ("leastsquares_helmholtz.e", context.equation_systems);
// 	}

// 	void mass_matrix_prod(LibMeshInit &init)
// 	{
// 		auto mesh = make_shared<libMesh::Mesh>(init.comm());		
// 		MeshTools::Generation::build_square (*mesh,
// 			10, 10,
// 			-1., 1.,
// 			-1., 1.,
// 			QUAD9);

// 		const int dim = mesh->mesh_dimension();

// 		LibMeshFEContext<LinearImplicitSystem> context(mesh);
		
// 		auto UXh = fe_space(LAGRANGE, FIRST, context);
// 		auto UYh = fe_space(LAGRANGE, FIRST, context);
		
// 		auto u_x = fe_function(UXh);
// 		auto u_y = fe_function(UYh);

// 		context.equation_systems.init();

// 		auto u  = prod(u_x, u_y);//, p);
// 		auto bf = integral(dot(u, u));

// 		u_x.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		u_y.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
		
// 		u_x.set_element(0);
// 		u_y.set_element(0);

// 		// p.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		// p.set_element(0);

// 		DenseMatrix<Real> mat;
// 		LibMeshBackend backend;
// 		backend.assemble(bf, mat);
// 		std::cout << "--------------------------------------------------------------------------\n" << std::endl;
// 		std::cout << mat.m() << " x " << mat.n() << "\n";
// 		mat.print();
// 		std::cout << "--------------------------------------------------------------------------\n" << std::endl;
// 	}

// 	void matrix_fe_tensor_product_fe_2x2(LibMeshInit &init)
// 	{
// 		auto mesh = make_shared<libMesh::Mesh>(init.comm());		
// 		MeshTools::Generation::build_square (*mesh,
// 			1, 1,
// 			0., 1.,
// 			0., 1.,
// 			TRI3);

// 		const int dim = mesh->mesh_dimension();

// 		auto order = FIRST;
// 		LibMeshFEContext<LinearImplicitSystem> context(mesh);
// 		auto UX = fe_space(LAGRANGE, order, context); auto u_x = fe_function(UX);
// 		auto UY = fe_space(LAGRANGE, order, context); auto u_y = fe_function(UY);
// 		auto UZ = fe_space(LAGRANGE, order, context); auto u_z = fe_function(UZ);
// 		auto UW = fe_space(LAGRANGE, order, context); auto u_w = fe_function(UW);

// 		context.equation_systems.init();

// 		u_x.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		u_y.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		u_z.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		u_w.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
		
// 		int index = 0;
// 		u_x.set_element(index);
// 		u_y.set_element(index);
// 		u_z.set_element(index);
// 		u_w.set_element(index);

// 		DenseVector<Real> b(dim);
// 		for(int i = 0; i < dim; ++i) { b(i) = 1.0; }

// 		auto u  = tensor_product_fe_2x2(u_x, u_y, u_z, u_w);
// 		auto f  = vec_coeff(b);

// 		const int dim2 = dim * dim;
// 		LMDenseMatrix A = identity(dim2, dim2);	
// 		// A.set(0, 0, 2);
// 		// A.set(dim2-1, dim2-1, 4);

// 		// A = make_stress_strain_rel_tensor(dim, 1, 1);


// 		auto bf = integral(dot(div(u), div(u)));
// 		auto mf = integral(dot(A * u, A * u));
// 		auto lf = integral(dot(f, div(u)));

// 		DenseMatrix<Real> mat, mat_m;
// 		DenseVector<Real> vec;
// 		LibMeshBackend backend;
// 		backend.assemble(bf, mat);
// 		backend.assemble(mf, mat_m);
// 		backend.assemble(lf, vec);


// 		std::cout << "--------------------------------------------------------------------------\n" << std::endl;
// 		std::cout << mat.m() << " x " << mat.n() << "\n";
// 		mat.print();
// 		std::cout << "--------------------------------------------------------------------------\n" << std::endl;
// 		mat_m.print();
// 		std::cout << "--------------------------------------------------------------------------\n" << std::endl;
// 		vec.print(std::cout);
// 		std::cout << "--------------------------------------------------------------------------\n" << std::endl;
// 	}

// 	void div_prod(LibMeshInit &init)
// 	{
// 		auto mesh = make_shared<libMesh::Mesh>(init.comm());		
// 		MeshTools::Generation::build_square (*mesh,
// 			10, 10,
// 			-1., 1.,
// 			-1., 1.,
// 			QUAD9);

// 		const int dim = mesh->mesh_dimension();

// 		LibMeshFEContext<LinearImplicitSystem> context(mesh);
// 		auto UX = fe_space(LAGRANGE, FIRST, context); 		auto u_x   = fe_function(UX);
// 		auto UY = fe_space(LAGRANGE, FIRST, context); 		auto u_y   = fe_function(UY);
// 		auto U  = vector_fe_space(LAGRANGE_VEC, FIRST, context); 	auto u_vec = fe_function(U);

// 		context.equation_systems.init();

// 		auto u  = prod(u_x, u_y);
// 		auto bf = integral(dot(div(u), div(u)));

// 		u_x.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		u_y.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		u_x.set_element(0);
// 		u_y.set_element(0);	

// 		u_vec.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));
// 		u_vec.set_element(0);

// 		DenseMatrix<Real> mat;
// 		LibMeshBackend backend;
// 		backend.assemble(bf, mat);

// 		std::cout << "--------------------------------------------------------------------------\n" << std::endl;
// 		mat.print();
// 		mat.resize(0, 0);
// 		std::cout << "--------------------------------------------------------------------------\n" << std::endl;
// 		auto bf_s = integral(dot(div(u_vec), div(u_vec)));
// 		backend.assemble(bf_s, mat);
// 		mat.print();
// 	}

// 	template<class U, class Sigma, class A11, class A12, class A21, class A22, class RHS1, class RHS2>
// 	void assemble_ls_elast(
// 		U &u,
// 		Sigma &sigma,
// 		const std::tuple<A11, A12, A21, A22> &blocks, 
// 		const std::tuple<RHS1, RHS2> &rhs_blocks, 
// 		LibMeshFEContext<LinearImplicitSystem> &context,
// 		const bool handle_boundary_conds) {

// 		const std::shared_ptr<MeshBase> &mesh = context.mesh;

// 		const auto &bf_11 = std::get<0>(blocks); 
// 		const auto &bf_12 = std::get<1>(blocks); 
// 		const auto &bf_21 = std::get<2>(blocks); 
// 		const auto &bf_22 = std::get<3>(blocks); 

// 		const auto &lf_1  = std::get<0>(rhs_blocks);
// 		const auto &lf_2  = std::get<1>(rhs_blocks);

// 		// Chrono c;
// 		auto ass = make_assembly([&]() -> void {

// 			// int i = 0;

// 			// c.start();
// 			DenseMatrix<Real> mat_11, mat_12, mat_21, mat_22;
// 			DenseVector<Real> vec_1, vec_2;
// 			std::vector<dof_id_type>  dof_u, dof_sigma;

// 			auto e_begin = mesh->active_local_elements_begin();
// 			auto e_end   = mesh->active_local_elements_end();

// 			LibMeshBackend backend;
// 			for(auto e_it = e_begin; e_it != e_end; ++e_it) {
// 				u.set_element(**e_it);

// 				for(int i = 0; i < sigma.size(); ++i) {
// 					sigma.get(i).set_element(**e_it);
// 				}

// 				backend.assemble(bf_11, mat_11);
// 				backend.assemble(bf_12, mat_12);
// 				backend.assemble(bf_21, mat_21);
// 				backend.assemble(bf_22, mat_22);

// 				backend.assemble(lf_1, vec_1);
// 				backend.assemble(lf_2, vec_2);

// 				// if(i == 0) {
// 				// 	std::cout << "------------------\n";
// 				// 	mat_11.print(std::cout);
// 				// 	std::cout << "------------------\n";
// 				// 	std::cout << "------------------\n";
// 				// 	mat_12.print(std::cout);
// 				// 	std::cout << "------------------\n";
// 				// 	std::cout << "------------------\n";
// 				// 	mat_21.print(std::cout);
// 				// 	std::cout << "------------------\n";
// 				// 	std::cout << "------------------\n";
// 				// 	mat_22.print(std::cout);
// 				// 	std::cout << "------------------\n";
// 				// 	i++;
// 				// }

// 				assert(is_symmetric(mat_11));
// 				assert(is_symmetric(mat_22));
// 				assert(is_approx_equal_tr(mat_21, mat_12));

// 				u.dof_map().dof_indices(*e_it, dof_u, u.var_num());

// 				if(handle_boundary_conds) {
// 					u.dof_map().constrain_element_vector(vec_1, dof_u);
// 				}

// 				//assemble in libmesh matrices
// 				if(handle_boundary_conds) {
// 					backend.constrain_mixed_matrix_and_add(*e_it, mat_11, mat_12, u, sigma, *context.system.matrix);
// 				} else {
// 					backend.mixed_matrix_and_add(*e_it, mat_11, mat_12, u, sigma, *context.system.matrix);
// 				}

// 				backend.mixed_matrix_add(*e_it, mat_21, mat_22, sigma, u, *context.system.matrix);
// 				backend.mixed_vector_add(*e_it, vec_1, vec_2, u, sigma, *context.system.rhs);
// 			}

// 			// c.stop();
// 			std::cout << "Assembly time: ";
// 			// c.describe(std::cout);
// 		});

// 		context.system.attach_assemble_object(ass);
// 		context.equation_systems.print_info();
// 		context.equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 0;
// 		context.equation_systems.solve();
// 	}

// 	template<class U, class Sigma>
// 	void assemble_ls_elast_system(U &u, Sigma &sigma, const double mu, const double lambda, 
// 								  DenseVector<Real> &force,
// 								  LibMeshFEContext<LinearImplicitSystem> &context, const bool handle_boundary_conds = true)
// 	{
// 		const int dim = context.mesh->mesh_dimension();

// 		LMDenseMatrix A = make_stress_strain_rel_tensor(dim, mu, lambda);

// 		DenseVector<Real> zero(dim);
// 		zero.zero();

// 		auto f0 = vec_coeff(zero);
// 		auto f  = vec_coeff(force);

// 		auto e    	   = 0.5 * (transpose(grad(u)) + grad(u)); 
		
// 		auto A_x_sigma = A * sigma;
// 		auto bf_11 	   = integral(dot(e, e));
// 		auto bf_12 	   = integral(-dot(A_x_sigma, e));
// 		auto bf_21 	   = integral(-dot(e, A_x_sigma));
// 		auto bf_22 	   = integral(dot(div(sigma), div(sigma)) + dot(A_x_sigma, A_x_sigma));

// 		auto lf_1      = integral(dot(f0, u));
// 		auto lf_2      = integral(-dot(f, div(sigma)));

// 		assemble_ls_elast(
// 			u, 
// 			sigma, 
// 			std::make_tuple(bf_11, bf_12, bf_21, bf_22), 
// 			std::make_tuple(lf_1, lf_2), context, handle_boundary_conds);
// 	}

// 	template<class U, class Sigma>
// 	void solve_ls_system(
// 		U &u, Sigma &sigma, 
// 		LibMeshFEContext<LinearImplicitSystem> &context,
// 		const std::string &output_path,
// 		const bool handle_boundary_conds = false) {

// 		const std::shared_ptr<MeshBase> &mesh = context.mesh;
// 		long n = context.system.solution->size();

// 		DSMatrixd mat;
// 		DVectord sol; 
// 		DVectord rhs; 
// 		mat = sparse(n, n, 40);
// 		sol = zeros(n);
// 		rhs = zeros(n);

// 		std::cout << "LibMesh solve terminated" << std::endl;
		
// 	#ifndef LIBMESH_HAVE_PETSC
// 		static_assert(false, "needs a libmesh-petsc installation");
// 	#endif	

// 		std::cout << "converting..." << std::endl;

// 		convert(*context.system.rhs, rhs);
// 		convert(*context.system.matrix, mat);

		

// 		if(handle_boundary_conds) {
// 			apply_boundary_conditions(u, mat, rhs);

// 			for(int i = 0; i < sigma.size(); ++i) {
// 				apply_boundary_conditions(sigma.get(i), mat, rhs);
// 			}

// 		}

// 		// write("mat.m", mat);

// 		std::cout << "starting solve" << std::endl;

// 		// Chrono c;
// 		// c.start();
		
// 		Factorization<DSMatrixd, DVectord>().solve(mat, rhs, sol);
// 		convert(sol, *context.system.solution);

// 		// disp(mat);
// 		// disp(rhs);
// 		// disp(sol);

// 		// c.stop();
// 		std::cout << "Solving time: ";
// 		// c.describe(std::cout);

// 		ExodusII_IO(*mesh).write_equation_systems (output_path, context.equation_systems);
// 	}

// 	void least_squares_elasticity_2D(LibMeshInit &init)
// 	{
// 		auto mesh = make_shared<libMesh::Mesh>(init.comm());		
// 		MeshTools::Generation::build_square (*mesh,
// 			10, 10,
// 			0, 1.,
// 			0, 1.,
// 			QUAD9);

// 		const int dim = mesh->mesh_dimension();

// 		LibMeshFEContext<LinearImplicitSystem> context(mesh);
		
// 		//stress
// 		auto Sigma = fe_tensor_product_space<2, 2>(LAGRANGE, FIRST, context);
// 		auto sigma = fe_function(Sigma);

// 		//displacement
// 		auto U = vector_fe_space("disp_", LAGRANGE_VEC, FIRST, context); 
// 		auto u = fe_function(U);

// 		std::function<void (const Point &, DenseVector<Real> &output)> left_bound_cond = [dim](const Point &, DenseVector<Real> &output) {
// 			output.resize(dim * dim + dim);
// 			output.zero();

// 			//offset to y coordinate of displacement
// 			output(dim * dim) = 0;
// 		};

// 		std::function<void (const Point &, DenseVector<Real> &output)> right_bound_cond = [dim](const Point &, DenseVector<Real> &output) {
// 			output.resize(dim * dim + dim);
// 			output.zero();

// 			//offset to y coordinate of displacement
// 			output(dim * dim) = 0.1;
// 			output(dim * dim + 1) = 0.2;
// 		};

// 		strong_enforce( boundary_conditions( u == vec_coeff(left_bound_cond),    { 1 }));
// 		strong_enforce( boundary_conditions( u == vec_coeff(right_bound_cond),   { 3 }));

// 		for(int i = 0; i < sigma.size(); ++i) {
// 			strong_enforce( boundary_conditions( sigma.get(i) == coeff(0.0), { 0, 2 }));
// 		}

// 		context.equation_systems.init();

// 		u.set_quad_rule(make_shared<libMesh::QGauss>(dim, SIXTH));

// 		for(int i = 0; i < sigma.size(); ++i) {
// 			sigma.get(i).set_quad_rule(make_shared<libMesh::QGauss>(dim, SIXTH));
// 		}
		
// 		double mu = 1.0, lambda = 1.0;

// 		DenseVector<Real> vec(dim);
// 		vec.zero();
// 		// vec(0) = -.8;
// 		// vec(1) = -.8;
 
// 		assemble_ls_elast_system(u, sigma, mu, lambda, vec, context, false);
// 		solve_ls_system(u, sigma, context,"least_squares_elasticity_2D.e", true);
// 	}

// 	void least_squares_elasticity_3D(LibMeshInit &init)
// 	{
// 		auto mesh = make_shared<libMesh::Mesh>(init.comm());		
// 		MeshTools::Generation::build_cube (*mesh,
// 			15, 15, 15,
// 			-1., 1., 
// 			-1., 1.,
// 			-1., 1.,
// 			HEX27);

// 		const int dim = mesh->mesh_dimension();

// 		LibMeshFEContext<LinearImplicitSystem> context(mesh);
		
// 		//stress
// 		auto Sigma = fe_tensor_product_space<3, 3>(LAGRANGE, FIRST, context);
// 		auto sigma = fe_function(Sigma);

// 		//displacement
// 		auto U = vector_fe_space("disp_", LAGRANGE_VEC, SECOND, context); 
// 		auto u = fe_function(U);

// 		std::function<void (const Point &, DenseVector<Real> &output)> cond_value = [&u](const Point &, DenseVector<Real> &output) {
// 			output.zero();
// 		};

// 		strong_enforce( boundary_conditions(u == vec_coeff(cond_value), {1, 3}) );

// 		for(int i = 0; i < sigma.size(); ++i) {
// 			strong_enforce( boundary_conditions( sigma.get(i) == coeff(0.0), { 0, 2, 4, 5 }));
// 		}


// 		context.equation_systems.init();

// 		u.set_quad_rule(make_shared<libMesh::QGauss>(dim, NINTH));

// 		for(int i = 0; i < sigma.size(); ++i) {
// 			sigma.get(i).set_quad_rule(make_shared<libMesh::QGauss>(dim, NINTH));
// 		}
		
// 		//bilinear forms
// 		double mu = 1.0, lambda = 1.0;

// 		DenseVector<Real> vec(dim);
// 		vec.zero();
// 		vec(1) = -.8;

// 		assemble_ls_elast_system(u, sigma, mu, lambda, vec, context, false);
// 		solve_ls_system(u, sigma, context,"least_squares_elasticity_3D.e", true);
// 	}

// 	void least_squares_contact(LibMeshInit &init)
// 	{
// 		std::cout << "-----------------------------\n";
// 		std::cout << "least_squares_contact\n";
// 		//////////////////////////////////////////////////
// 		//////////////////////////////////////////////////

// 		ContactSimParams params = contact_least_squares_2;

// 		auto mesh = make_shared<libMesh::Mesh>(init.comm());		
// 		mesh->read(params.mesh_path);
// 		plot_mesh(*mesh, "mesh");

// 		const int dim = mesh->mesh_dimension();

// 		assert(dim == 2);

// 		LibMeshFEContext<LinearImplicitSystem> context(mesh);

// 		//stress
// 		auto Sigma = fe_tensor_product_space<2, 2>(LAGRANGE, FIRST, context);
// 		auto sigma = fe_function(Sigma);

// 		//displacement
// 		auto U = vector_fe_space("disp_", LAGRANGE_VEC, FIRST, context); 
// 		auto u = fe_function(U);


// 		std::function<void (const Point &, DenseVector<Real> &output)> top_bound_cond = [dim,&params](const Point &, DenseVector<Real> &output) {
// 			output.resize(dim * dim + dim);
// 			output.zero();

// 			//offset to y coordinate of displacement
// 			output(dim * dim + 1) = params.dirichlet_value_1;
// 		};

// 		std::function<void (const Point &, DenseVector<Real> &output)> bottom_bound_cond = [dim,&params](const Point &, DenseVector<Real> &output) {
// 			output.resize(dim * dim + dim);
// 			output.zero();

// 			//offset to y coordinate of displacement
// 			output(dim * dim + 1) = params.dirichlet_value_2;
// 		};
		
// 		auto bc_1 = boundary_conditions( u == vec_coeff(top_bound_cond),    { params.boundary_tag_1 } );
// 		auto bc_2 = boundary_conditions( u == vec_coeff(bottom_bound_cond), { params.boundary_tag_2 } );

// 		strong_enforce(bc_1);
// 		strong_enforce(bc_2);

// 		for(int i = 0; i < sigma.size(); ++i) {
// 			strong_enforce( boundary_conditions( sigma.get(i) == coeff(0.0), { 3 }));
// 		}

// 		context.equation_systems.init();

// 		u.set_quad_rule(make_shared<libMesh::QGauss>(dim, SIXTH));


// 		for(int i = 0; i < sigma.size(); ++i) {
// 			sigma.get(i).set_quad_rule(make_shared<libMesh::QGauss>(dim, SIXTH));
// 		}
		
// 		double mu = 1.0, lambda = 1.0;

// 		DenseVector<Real> vec(dim);
// 		vec.zero();

// 		assemble_ls_elast_system(u, sigma, mu, lambda, vec, context, false);

// 		//////////////////////////////////////////////////
// 		//////////////////////////////////////////////////

// 		// Chrono c;
// 		// c.start();

// 		long n = u.dof_map().n_dofs();

// 		MortarContactAssembler assembler(make_ref(U));

// 		DSMatrixd coupling;
// 		DVectord gap;
// 		DVectord normals;
// 		DSMatrixd orhtogonal_trafos;

// 		std::vector<bool> is_contact_node;
// 		if(!assembler.assemble(coupling, 
// 			gap, 
// 			normals, 
// 			orhtogonal_trafos, 
// 			is_contact_node, 
// 			params.search_radius))
// 		{
// 			//Just set some values that do not change the original system
// 			coupling = identity(n, n);
// 			orhtogonal_trafos = identity(n, n);
// 			gap = values(n, 10000);
// 			normals = zeros(n);
// 			is_contact_node.resize(n);
// 			std::fill(is_contact_node.begin(), is_contact_node.end(), false);
// 		}

// 		DSMatrixd K  = sparse(n, n, 20);
// 		DVectord  rhs = zeros(n);

// 		utopia::convert(*context.system.rhs, rhs);
// 		utopia::convert(*context.system.matrix, K);


// 		disp(size(coupling));
// 		disp(size(orhtogonal_trafos));
// 		disp(size(K));
// 		disp(is_contact_node.size());

// 		// disp(gap);
// 		// disp(normals);
// 		// print_vector(is_contact_node.begin(), is_contact_node.end());
// 		// disp((coupling));
// 		// disp((orhtogonal_trafos));
// 		// disp((K));

// 		//Change of basis
// 		DVectord  sol_c = zeros(size(rhs));
// 		DVectord  rhs_c = transpose(orhtogonal_trafos) * transpose(coupling) * rhs;
// 		DSMatrixd K_c   = transpose(orhtogonal_trafos) * DSMatrixd(transpose(coupling) * K * coupling) * orhtogonal_trafos;
// 		DVectord  gap_c = transpose(coupling) * gap;
// 		apply_boundary_conditions(u, K_c, rhs_c);

// 		// if(handle_boundary_conds) {
// 		// 	apply_boundary_conditions(u, mat, rhs);

// 		// 	for(int i = 0; i < sigma.size(); ++i) {
// 		// 		apply_boundary_conditions(sigma.get(i), mat, rhs);
// 		// 	}

// 		// }

// 		SemismoothNewton<DSMatrixd, DVectord> newton(std::make_shared<Factorization<DSMatrixd, DVectord> >());
// 		// SemismoothNewton<DSMatrixd, DVectord> newton(std::make_shared<ConjugateGradient<DSMatrixd, DVectord> >());

// 		newton.verbose(true);
// 		//newton.solve(sol_c, K_c, rhs_c, gap_c);

// 		//Change back to original basis
// 		DVectord sol = coupling * (orhtogonal_trafos * sol_c);
		
// 		//plot stuff
// 		convert(sol, *context.system.solution);
// 		ExodusII_IO(*mesh).write_equation_systems ("least_squares_contact.e", context.equation_systems);

// 		convert(gap, *context.system.solution);
// 		ExodusII_IO(*mesh).write_equation_systems ("least_squares_contact_gap.e", context.equation_systems);

// 		convert(normals, *context.system.solution);
// 		ExodusII_IO(*mesh).write_equation_systems ("least_squares_contact_normals.e", context.equation_systems);
// 	}

// 	void run_least_squares_examples(libMesh::LibMeshInit  &init)
// 	{
// 		// mass_matrix_prod(init);
// 		// div_prod(init);
// 		// matrix_fe_tensor_product_fe_2x2(init);
// 		// leastsquares_helmoholtz(init);
// 		// least_squares_elasticity_2D(init);
// 		// least_squares_elasticity_3D(init);
// 		least_squares_contact(init);
// 	}
// }

#include <iostream>
#include <cmath>
#include "utopia_FEDSLBaseExamples.hpp"
#include "utopia.hpp"

//fe extension
#include "utopia_fe_core.hpp"
#include "MortarAssembler.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>

#include "utopia_UGMeshReader.hpp"


using namespace std;
using namespace libMesh;

namespace utopia {

	static void upper_boundary_cond(const Point & p, DenseVector<Real> & output)
	{
		output.resize(3);
		output.zero();
		output(1) = -0.1;
	}

	static void lower_boundary_cond(const Point & p, DenseVector<Real> & output)
	{
		output.resize(3);
		output.zero();
		output(1) = 0.1;
	}

	void anisotropic_laplacian(LibMeshInit &init)
	{
		auto mesh = make_shared<Mesh>(init.comm());		
		MeshTools::Generation::build_square (*mesh,
			15, 15,
			-1., 1.,
			-1., 1.,
			QUAD9);

		for(auto it = mesh->active_local_elements_begin(); it != mesh->active_local_elements_end(); ++it) {
			auto &n = (*it)->node_ref(0);

			if(n(0) < -0.5) {
				(*it)->subdomain_id() = 1;
			} else if(n(0) < 0.5) {
				(*it)->subdomain_id() = 2;
			} else {
				(*it)->subdomain_id() = 3;
			}
		}

		LibMeshFEContext<LinearImplicitSystem> context(mesh);
		auto Vh = fe_space(LAGRANGE, FIRST, context);
		auto v  = fe_function(Vh);

		strong_enforce( boundary_conditions(v == coeff(1), {2}) );
		context.equation_systems.init();

		const int dim = mesh->mesh_dimension();
		v.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));

		//utopia Wrapper implemented with inliner
		LMDenseMatrix A1 = identity(dim, dim);
		{
			Write<LMDenseMatrix> w_A(A1);
			A1.set(0, 0, 10);
			A1.set(1, 1, 0.1);
			if(dim > 2) A1.set(2, 2, 1);
		}

		LMDenseMatrix A2 = identity(dim, dim);
		{
			Write<LMDenseMatrix> w_A(A2);
			A2.set(0, 0, 0.1);
			A2.set(1, 1, 10);
			if(dim > 2) A2.set(2, 2, 1);
		}

		LMDenseMatrix A_default = identity(dim, dim);
		{
			Write<LMDenseMatrix> w_A(A_default);
			A_default.set(0, 0, 1);
			A_default.set(1, 1, 1);
			if(dim > 2) A_default.set(2, 2, 1);
		}

		auto A = block_var(A_default, {{1, A1}, {2, A2}});

		std::function<Real (const Point &)> f = [dim](const Point &p) {
			return 10 * std::sqrt( p(0) * p(0) + p(1) * p(1) ) - 5.0;
		};

		auto ass = make_assembly([&]() -> void {
			assemble(v, v, integral(dot(A * grad(v), grad(v)) ), integral(dot(coeff(f), v)), 
				*context.system.matrix,  
				*context.system.rhs);
		});

		context.system.attach_assemble_object(ass);
		// context.equation_systems.print_info();
		context.equation_systems.solve();

		DVectord u_sol;
		convert(*context.system.solution, u_sol);



		ExodusII_IO(*mesh).write_equation_systems ("anisotropic_laplacian.e", context.equation_systems);

		double *arr;
		VecGetArray(raw_type(u_sol), &arr);
		plot_mesh_f(v.mesh(), arr, "anisotropic_laplacian");
		VecRestoreArray(raw_type(u_sol), &arr);
	}

	void linear_elasticity(LibMeshInit &init)
	{

		auto mesh = make_shared<Mesh>(init.comm());		
		MeshTools::Generation::build_cube (*mesh,
			10, 10, 10,
			-1., 1.,
			-1., 1.,
			-1., 1.,
			HEX27);

		const int dim = mesh->mesh_dimension();

		LibMeshFEContext<LinearImplicitSystem> context(mesh);
		auto Uh = vector_fe_space("disp_", LAGRANGE_VEC, SECOND, context);
		auto u  = fe_function(Uh);	

		strong_enforce( boundary_conditions(u == vec_coeff(upper_boundary_cond), {3}) );
		strong_enforce( boundary_conditions(u == vec_coeff(lower_boundary_cond), {1}) );

	//system initialized
		context.equation_systems.init();

	//FIXME to be moved in the assembly loop
		u.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));

	//bilinear forms
	double mu = 1.0, lambda = 1.0;
	auto e      = transpose(grad(u)) + grad(u); //0.5 moved below -> (2 * 0.5 * 0.5 = 0.5)
	auto b_form = integral((mu * 0.5) * dot(e, e) + lambda * dot(div(u), div(u)));

	DenseVector<Real> vec(dim);
	vec.zero();
	vec(2) = -0.8;

	//linear forms
	auto f 	    = vec_coeff(vec);
	auto l_form = integral(dot(f, u));

	auto &matrix = *context.system.matrix;
	auto &rhs 	 = *context.system.rhs;

	auto ass = make_assembly([&]() -> void {
		double t = MPI_Wtime();
		
		assemble(u, u, b_form, l_form, matrix, rhs);

		t = MPI_Wtime() - t;

		printf("--------------------------------\n");
		printf("Assembly: %g seconds\n", t);
		printf("--------------------------------\n");
	});



	context.system.attach_assemble_object(ass);
	// context.equation_systems.print_info();
	context.equation_systems.solve();

	ExodusII_IO(*mesh).write_equation_systems ("linear_elasticity.e", context.equation_systems);

}

void nonlinear_laplace_eq(LibMeshInit &init)
{
	int n = 15;
	auto mesh = make_shared<Mesh>(init.comm());		
	MeshTools::Generation::build_square (*mesh,
		n, n,
		-1., 1.,
		-1., 1.,
		QUAD9);


	for(auto it = mesh->active_local_elements_begin(); it != mesh->active_local_elements_end(); ++it) {
		auto &n = (*it)->node_ref(0);

		// std::cout << (*it)->subdomain_id() << std::endl;

		if(n(0) < 0.0) {
			(*it)->subdomain_id() = 1;
		} else {
			(*it)->subdomain_id() = 2;
		}
	}

	LibMeshFEContext<LinearImplicitSystem> context(mesh);
	auto Vh = fe_space(LAGRANGE, FIRST, context);
	auto v  = fe_function(Vh);

	strong_enforce( boundary_conditions(v == coeff(-0.1), {0, 1, 2, 3}) );

	context.equation_systems.init();

	const int dim = mesh->mesh_dimension();
	v.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));

	std::function<Real(const Point &p)> rhs_fun_1 = [dim](const Point &p) -> Real {
		Real center = 0;

		Real ret = 0;
		for(int i = 0; i < dim; ++i) {
			ret += (p(i) - center) * (p(i) - center);
		}

		return 10 * std::sqrt(ret) - 5;
	};

	std::function<Real(const Point &p)> rhs_fun_2 = [dim](const Point &) -> Real {
		return 1;
	};

	auto c1 = coeff(rhs_fun_1);
	auto c2 = coeff(rhs_fun_2);

	LMDenseMatrix A = identity(dim, dim);
	{	
		Write<LMDenseMatrix> w_A(A);
		A.set(0, 0, 2.);
		A.set(1, 1, 0.5);
		if(dim > 2) A.set(2, 2, 1);
	}

	auto v_k = interpolate(coeff(1.0), v, make_ref(*context.system.solution));
	auto bf  = integral(dot(1./(pow2(v_k) + coeff(0.1)) * (A * grad(v)), grad(v)));	
	auto lf  = integral(dot(c1, v), 1) + integral(dot(c2, v), 2);

	auto ass = make_assembly([&]() -> void {
		std::cout << "assemble called, norm(u):" << context.system.solution->l2_norm() << std::endl;
		assemble(v, v, bf, lf, *context.system.matrix,  *context.system.rhs);
	});

	context.system.attach_assemble_object(ass);
	// context.equation_systems.print_info();

	for(auto i = 0; i < 10; ++i) {
		context.equation_systems.solve();
	}

	ExodusII_IO(*mesh).write_equation_systems ("nonlinear_laplace_eq.e", context.equation_systems);
}

void boundary_conds(LibMeshInit &init)
{
	auto mesh = make_shared<Mesh>(init.comm());		
	MeshTools::Generation::build_square (*mesh,
		10, 10,
		-1., 1.,
		-1., 1.,
		QUAD9);

	const int dim = mesh->mesh_dimension();

	LibMeshFEContext<LinearImplicitSystem> context(mesh);

	auto Uh = fe_space(LAGRANGE, FIRST, context); 
	auto u  = fe_function(Uh);


	std::function<Real(const Point &)> rhs_fun = [dim](const Point &p) -> Real {
		Real center = 0;
		Real ret = 0;
		for(int i = 0; i < dim; ++i) {
			ret += (p(i) - center) * (p(i) - center);
		}

		ret = 10 * std::sqrt(ret) - 5;
		return ret;
	};

	auto g  = coeff(rhs_fun);

	strong_enforce( boundary_conditions(u == g, {0, 1}) );

	//examples for specifying variational inequalities:
	//bc = boundary_conditions(normal_projection(u, {0})  <= gap_to(surf_mesh, {1}));
	//bc = boundary_conditions(normal_projection(u, {0})  <= gap([](Point &origin, Point &direction) -> double { /*return intersection*/}));
	//c  = boundary_constraint(normal_projection(jump(u)) <= gap()));
	//c  = constraint(lower <= u <= upper);
}

void run_base_examples(libMesh::LibMeshInit &init)
{
	anisotropic_laplacian(init);
	linear_elasticity(init);
	nonlinear_laplace_eq(init);
	boundary_conds(init);
}

}
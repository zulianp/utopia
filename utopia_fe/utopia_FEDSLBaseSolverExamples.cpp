
#include <iostream>
#include <cmath>
#include "utopia_FEDSLBaseSolverExamples.hpp"
#include "utopia.hpp"

//fe extension
#include "utopia_fe.hpp"
// #include "MortarAssembler.hpp"
#include "utopia_Socket.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>


#include "libmesh/nonlinear_implicit_system.h"


using namespace std;
using namespace libMesh;

namespace utopia {

	// static void upper_boundary_cond(const Point & p, DenseVector<Real> & output)
	// {
	// 	output.resize(3);
	// 	output.zero();
	// 	output(1) = -0.1;
	// }

	// static void lower_boundary_cond(const Point & p, DenseVector<Real> & output)
	// {
	// 	output.resize(3);
	// 	output.zero();
	// 	output(1) = 0.1;
	// }



// void nonlinear_laplace_eq_solver_test(LibMeshInit &init)
// {
// 	int n = 50;
// 	auto mesh = make_shared<Mesh>(init.comm());		
// 	MeshTools::Generation::build_square (*mesh,
// 		n, n,
// 		-1., 1.,
// 		-1., 1.,
// 		QUAD9);

	
// 	LibMeshFEContext<LinearImplicitSystem> context(mesh);
// 	// LibMeshFEContext<NonlinearImplicitSystem> context(mesh); 
// 	auto Vh = fe_space(LAGRANGE, FIRST, context);
// 	auto v  = fe_function(Vh);

// 	strong_enforce( boundary_conditions(v == coeff(-0.1), {0, 1, 2, 3}) );

// 	context.equation_systems.init();

// 	const int dim = mesh->mesh_dimension();
// 	v.set_quad_rule(make_shared<libMesh::QGauss>(dim, FIFTH));

// 	std::function<Real(const Point &p) > rhs_fun = [dim](const Point &p) -> Real {
// 		Real center = 0;

// 		Real ret = 0;
// 		for(int i = 0; i < dim; ++i) {
// 			ret += (p(i) - center) * (p(i) - center);
// 		}

// 		return 10 * std::sqrt(ret) - 5;
// 	};

// 	auto c  = coeff(rhs_fun);
	
// 	LMDenseMatrix A = identity(dim, dim);	
// 	A.set(0, 0, 2.);
// 	A.set(1, 1, 0.5);
// 	if(dim > 2) A.set(2, 2, 1);

// 	auto v_k = interpolate(coeff(1.0), v, make_ref(*context.system.solution));
// 	auto bf  = integral(dot(1./(pow2(v_k) + coeff(0.1)) * (A * grad(v)), grad(v)));	
// 	auto lf  = integral(dot(c, v));

// 	auto ass = make_assembly([&]() -> void {
// 		// std::cout << "assemble called, norm(u):" << context.system.solution->l2_norm() << std::endl;
// 		assemble(v, v, bf, lf, *context.system.matrix,  *context.system.rhs);
// 	});

// 	context.system.attach_assemble_object(ass);
// 	context.equation_systems.print_info();

// 	// for(auto i = 0; i < 10; ++i) {
// 		context.equation_systems.solve();
// 	// }

// 	ExodusII_IO(*mesh).write_equation_systems ("nonlinear_laplace_eq.e", context.equation_systems);
// }



	// void run_solver_ex(libMesh::LibMeshInit &init)
	// {

	// 	nonlinear_laplace_eq_solver_test(init);
	// }

}
#include "utopia.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/linear_partitioner.h>

//fe extension
#include "utopia_fe.hpp"
#include "MortarAssembler.hpp"
#include "MixedParMortarAssembler.hpp"
#include "ParMortarAssembler.hpp"
#include "utopia_Socket.hpp"
#include "utopia_ContactSimParams.hpp"

using namespace utopia;
using namespace std;
using namespace libMesh;


template<class FE>
void assemble_linear_elasticity(FE &fe, const Real lambda, const Real mu, libMesh::DenseMatrix<libMesh::Real> &mat)
{
	typedef libMesh::DenseMatrix<libMesh::Real> DenseMatrixT;
	typedef libMesh::TensorValue<libMesh::Real> TensorValueT;
	typedef libMesh::DenseVector<libMesh::Real> DenseVectorT;
	typedef unsigned int uint;

	auto &lm_fe = fe.get_fe();
	// auto &fun 	= lm_fe.get_phi();
	auto &grad	= lm_fe.get_dphi();
	auto &JxW	= lm_fe.get_JxW();
	auto &div   = lm_fe.get_div_phi();

	uint n_quad_points = grad[0].size();

	mat.resize(grad.size(), grad.size());
	mat.zero();

	std::vector<TensorValueT> strain(grad.size());

	for (uint qp = 0; qp < n_quad_points; qp++) {
		//precompute straint tensor for each quadrature point
		for(uint i = 0; i < grad.size(); ++i) {
			strain[i]  = grad[i][qp];
			strain[i] += grad[i][qp].transpose();
		}

		for (uint i = 0; i < strain.size(); i++) {
			for (uint j = i; j < strain.size(); j++) {
				mat(i, j) += ( mu * 0.5 * strain[i].contract(strain[j]) + lambda * div[i][qp] * div[j][qp] ) * JxW[qp];
			}
		}
	}

	//exploit symmetry
	for(uint i = 0; i < mat.n(); ++i) {
		for(uint j = i+1; j < mat.n(); ++j) {
			mat(j, i) = mat(i, j);
		}
	}
}

template<class FE>
void assemble_rhs(FE &fe, libMesh::DenseVector<libMesh::Real> &vec)
{
	typedef libMesh::DenseMatrix<libMesh::Real> DenseMatrixT;
	typedef libMesh::TensorValue<libMesh::Real> TensorValueT;
	typedef libMesh::DenseVector<libMesh::Real> DenseVectorT;
	typedef unsigned int uint;

	auto &lm_fe = fe.get_fe();
	auto &fun 	= lm_fe.get_phi();

	vec.resize(fun.size());
	vec.zero();
}

template<class FE, class Matrix, class Vector>
void assemble_linear_elasticity_system(
	const Real lambda, const Real mu, 
	FE &fe, 
	Matrix &mat,
	Vector &rhs)
{
	using namespace libMesh;
	LibMeshBackend backend;

	auto e_begin = fe.mesh().active_local_elements_begin();
	auto e_end   = fe.mesh().active_local_elements_end();

	std::vector<dof_id_type> dof_indices;

	DenseMatrix<Real> el_mat;
	DenseVector<Real> el_vec;
	for(auto e_it = e_begin; e_it != e_end; ++e_it) {
		fe.set_element(**e_it);


		assemble_linear_elasticity(fe, lambda, mu, el_mat);
		assemble_rhs(fe, el_vec);

		fe.dof_map().dof_indices(*e_it, dof_indices, fe.var_num());
		fe.dof_map().heterogenously_constrain_element_matrix_and_vector(el_mat, el_vec, dof_indices);			
		
		add_matrix(el_mat, dof_indices, dof_indices, mat);
		add_vector(el_vec, dof_indices, rhs);
	}
}



void run_biomechanics_example(libMesh::LibMeshInit &init)
{
	std::cout << "-----------------------------\n";
	std::cout << "run_biomechanics_example\n";
			//////////////////////////////////////////////////
			//////////////////////////////////////////////////

			// static const bool is_leaflet = true;
			// ContactSimParams params = leaflets_contact;

	static const bool is_leaflet = false;
			 // ContactSimParams params = contact_cuboids;
	    	// ContactSimParams params = contact8;
			// ContactSimParams params = triple_contact_circle;
			ContactSimParams params = multi_contact_3D_2;
	// ContactSimParams params = contact_cylinder;


	auto mesh = make_shared<Mesh>(init.comm());		
	mesh->read(params.mesh_path);
	plot_mesh(*mesh, "mesh");

	const int dim = mesh->mesh_dimension();

	//////////////////////////////////////////////////
	//////////////////////////////////////////////////

	Chrono c;
	c.start();

	LibMeshFEContext<LinearImplicitSystem> context(mesh);
	auto space = vector_fe_space("disp_", LAGRANGE_VEC, FIRST, context);
	auto u = fe_function(space);

	//boundary conditions
	std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_1 = [dim, &params](const Point &, DenseVector<Real> &output) {
		output.resize(dim);
		output.zero();

		if(is_leaflet) {
			// output(0) = -params.dirichlet_value_1;
			output(1) = params.dirichlet_value_1;
		} else {
			//offset to y coordinate of displacement
			output(1) = params.dirichlet_value_1;
		}
	};

	auto bc_1 = boundary_conditions(u == vec_coeff(boundary_cond_1),    { params.boundary_tag_1 });
	strong_enforce(bc_1);

	std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_2 = [dim, &params](const Point &, DenseVector<Real> &output) {
		output.resize(dim);
		output.zero();

		if(is_leaflet) {
			output(0) = -params.dirichlet_value_2;
			output(1) = -params.dirichlet_value_2;
		} else {
			//offset to y coordinate of displacement
			output(1) = params.dirichlet_value_2;
		}
	};

	auto bc_2 = boundary_conditions(u == vec_coeff(boundary_cond_2), { params.boundary_tag_2 });
	strong_enforce(bc_2);

	std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_3 = [dim, &params](const Point &, DenseVector<Real> &output) {
		output.resize(dim);
		output.zero();

		if(is_leaflet) {
			output(0) = params.dirichlet_value_3;
			output(1) = -params.dirichlet_value_3;
		} else {
			//offset to y coordinate of displacement
			output(1) = params.dirichlet_value_3;
		}
	};

	if(params.boundary_tag_3 >= 0) {
		auto bc_3 = boundary_conditions(u == vec_coeff(boundary_cond_3), { params.boundary_tag_3 });
		strong_enforce(bc_3);
	}

	context.equation_systems.init();
	u.set_quad_rule(make_shared<libMesh::QGauss>(dim, SIXTH));

	double mu   = 1.0, lambda = 1.0;
	DenseVector<Real> vec(dim);
	vec.zero();

	auto ass = make_assembly([&]() -> void {
		assemble_linear_elasticity_system(mu, lambda, u, *context.system.matrix, *context.system.rhs);
	});

	context.system.attach_assemble_object(ass);
	context.equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1;
	context.equation_systems.solve();

	c.stop();
	std::cout << "Matrix and rhs assembly: ";
	c.describe(std::cout); 
	std::cout << std::endl;

	long n = u.dof_map().n_dofs();

	MortarContactAssembler assembler(make_ref(space));
	assembler.set_strict_gap_policy(params.strict_search_radius);

	DSMatrixd coupling;
	DVectord gap;
	DVectord normals;
	DSMatrixd orhtogonal_trafos;

	std::vector<bool> is_contact_node;
	if(!assembler.assemble(coupling, 
		gap, 
		normals, 
		orhtogonal_trafos, 
		is_contact_node, 
		params.search_radius))
	{
		//Just set some values that do not change the original system
		coupling = identity(n, n);
		orhtogonal_trafos = identity(n, n);
		gap = values(n, 10000);
		normals = zeros(n);
		is_contact_node.resize(n);
		std::fill(is_contact_node.begin(), is_contact_node.end(), false);
	}

	DSMatrixd K  = sparse(n, n, 20);
	DVectord rhs = zeros(n);

	convert(*context.system.rhs, rhs);
	convert(*context.system.matrix, K);

	//Change of basis
	DVectord  sol_c = zeros(size(rhs));
	DVectord  rhs_c = transpose(orhtogonal_trafos) * transpose(coupling) * rhs;
	DSMatrixd K_c   = transpose(orhtogonal_trafos) * DSMatrixd(transpose(coupling) * K * coupling) * orhtogonal_trafos;
	
	// DVectord  gap_c = transpose(coupling) * gap;
	DVectord  gap_c = gap;
	apply_boundary_conditions(u, K_c, rhs_c);

	SemismoothNewton<DSMatrixd, DVectord> newton(std::make_shared<Factorization<DSMatrixd, DVectord> >());
	newton.verbose(true);
	newton.solve(sol_c, K_c, rhs_c, gap_c);

	//Change back to original basis
	DVectord sol = coupling * (orhtogonal_trafos * sol_c);

	convert(sol, *context.system.solution);
	ExodusII_IO(*mesh).write_equation_systems ("elasticity_contact.e", context.equation_systems);

	{
		Read<DVectord> w_g(gap);
		Read<DVectord> w_n(normals);

		for(size_t i = 0; i < is_contact_node.size(); ++i) {
			if(is_contact_node[i]) {
				const double len = gap.get((i/dim)*dim);
				context.system.solution->set(i, len * normals.get(i));
			} else {
				context.system.solution->set(i, 0);
			}
		}
	}

	ExodusII_IO(*mesh).write_equation_systems ("elasticity_contact_cn.e", context.equation_systems);
}
#include "utopia.hpp"
#include "utopia_assemble_contact.hpp"

#include <libmesh/const_function.h>
#include <libmesh/petsc_vector.h>
#include <libmesh/petsc_matrix.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/linear_partitioner.h>
#include <libmesh/mesh_refinement.h>

//fe extension
#include "utopia_fe.hpp"
#include "MortarAssembler.hpp"
#include "utopia_Socket.hpp"
#include "utopia_ContactSimParams.hpp"
#include "utopia_UGMeshReader.hpp"
#include "express_Communicator.hpp"




using namespace utopia;
using namespace std;
using namespace libMesh;

class LameeParameters {
public:
	LameeParameters(const double default_mu, const double default_lambda)
	: default_mu_(default_mu), default_lambda_(default_lambda)
	{}
	
	double mu(const int id) const
	{
		auto it = mu_.find(id);
		if(it == mu_.end()) {
			return default_mu_;
		}
		
		return it->second;
	}
	
	double lambda(const int id) const
	{
		auto it = lambda_.find(id);
		if(it == lambda_.end()) {
			return default_lambda_;
		}
		
		return it->second;
	}
	
	double default_mu_, default_lambda_;
	std::map<int, double> mu_;
	std::map<int, double> lambda_;
};

double von_mises_stress_2(const double *stress)
{
	double result =  0.5 * ( stress[0] - stress[3] ) *
	( stress[0] - stress[3] ) +
	3.0  *  stress[1] * stress[1];
	
	result = sqrt( fabs(result) );
	assert(result == result && "von_mises_stress_2: result is nan");
	return result;
}

double von_mises_stress_3(const double *stress)
{
	double result =  0.5 * ( stress[0] - stress[4] ) *
	( stress[0] - stress[4] ) +
	3.0  *  stress[1] * stress[1];
	
	result += 0.5 * (stress[8] - stress[4]) * (stress[8] - stress[4]) + 3.0  * stress[7] * stress[7];
	result += 0.5 * (stress[8] - stress[0]) * (stress[8] - stress[0]) + 3.0  * stress[6] * stress[6];
	
	result = sqrt( fabs(result) );
	
	assert(result == result && "von_mises_stress_3: result is nan");
	return result;
}

double von_mises_stress(const int n_dims, const double * stress)
{
	switch(n_dims) {
		case 2: { return von_mises_stress_2(stress); }
		case 3: { return von_mises_stress_3(stress); }
		default : { assert(false && "von_mises_stress: not supported for dim."); return 0.0; }
	}
	return 0.0;
}


void stress_linear_elasticity(const double mu, const double lambda, const DenseMatrix<Real> &grad_u, DenseMatrix<Real> &stress)
{
	stress = grad_u;
	const int n = stress.m();
	
	double trace_grad_u = 0;
	
	for(int i = 0; i < n; ++i) {
		trace_grad_u += grad_u(i, i);
		
		for(int j = 0; j < n; ++j) {
			stress(i, j) += grad_u(j, i);
		}
	}
	
	stress *= mu;
	double temp = (lambda * trace_grad_u);
	for(int i = 0; i < n; ++i) {
		stress(i, i) += temp;
	}
}

template<class FE>
void von_mises_stress_linear_elasticity(FE &fe, const int dims, const double mu, const double lambda, const libMesh::DenseVector<libMesh::Real> &u,
										libMesh::DenseVector<libMesh::Real> &von_mises_stress_vec,
										libMesh::DenseVector<libMesh::Real> &mass_vec)
{
	auto &fun = fe.get_fe().get_phi();
	auto &grad_fun = fe.get_fe().get_dphi();
	auto &JxW	= fe.get_fe().get_JxW();
	
	
	
	von_mises_stress_vec.resize(fun.size());
	von_mises_stress_vec.zero();
	
	mass_vec.resize(fun.size());
	mass_vec.zero();
	
	Real vm_stress = 0.0;
	Real mass = 0.0;
	
	DenseMatrix<Real> grad_u;
	DenseMatrix<Real> stress;
	for(SizeType qp = 0; qp < JxW.size(); ++qp) {
		DenseMatrix<Real> grad_u(dims, dims);
		grad_u.zero();
		
		for(SizeType i = 0; i < grad_fun.size(); ++i) {
			for(SizeType di = 0; di < dims;  ++di) {
				for(SizeType dj = 0; dj < dims; ++dj) {
					grad_u(di, dj) += grad_fun[i][qp](di, dj) * u(i);
				}
			}
			
		}
		
		stress_linear_elasticity(mu, lambda, grad_u, stress);
		//von mises
		
		for(SizeType i = 0; i < fun.size(); ++i) {
			auto FxJxW = fun[i][qp] * von_mises_stress(dims, &stress.get_values()[0]) * JxW[qp];
			auto MxJxW = fun[i][qp] * JxW[qp];
			
			for(SizeType d = 0; d < dims; ++d) {
				von_mises_stress_vec(i) += FxJxW(d);
				mass_vec(i) += MxJxW(d);
			}
		}
	}
}



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
									   const LameeParameters &params,
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
		
		const int block_id = (*e_it)->subdomain_id();
		
		const double mu 	= params.mu(block_id);
		const double lambda = params.lambda(block_id);
		
		assemble_linear_elasticity(fe, lambda, mu, el_mat);
		assemble_rhs(fe, el_vec);
		
		fe.dof_map().dof_indices(*e_it, dof_indices, fe.var_num());
		fe.dof_map().heterogenously_constrain_element_matrix_and_vector(el_mat, el_vec, dof_indices);
		
		add_matrix(el_mat, dof_indices, dof_indices, mat);
		add_vector(el_vec, dof_indices, rhs);
	}
}

template<class FE, class Vector>
void assemble_von_mises_stress(
							   const LameeParameters &params,
							   FE &fe,
							   const Vector &u,
							   Vector &stress)
{
	using namespace libMesh;
	LibMeshBackend backend;
	
	auto e_begin = fe.mesh().active_local_elements_begin();
	auto e_end   = fe.mesh().active_local_elements_end();
	
	std::vector<dof_id_type> dof_indices;
	
	Vector mass = local_zeros(local_size(u));
	stress = local_zeros(local_size(u));
	
	Read<Vector> r_u(u);
	{
		Write<Vector> w_s(stress), w_m(mass);
		
		DenseVector<Real> u_local, stress_local, mass_local;
		for(auto e_it = e_begin; e_it != e_end; ++e_it) {
			fe.set_element(**e_it);
			
			const int block_id = (*e_it)->subdomain_id();
			
			const double mu 	= params.mu(block_id);
			const double lambda = params.lambda(block_id);
			
			fe.dof_map().dof_indices(*e_it, dof_indices, fe.var_num());
			
			get_vector(u, dof_indices, u_local);
			
			von_mises_stress_linear_elasticity(fe, fe.mesh().mesh_dimension(), mu, lambda, u_local, stress_local, mass_local);
			add_vector(stress_local, dof_indices, stress);
			add_vector(mass_local, dof_indices, mass);
		}
	}
	
	stress = e_mul(stress, 1./mass);
}


static void plot_scaled_normal_field(MeshBase &mesh,
									 const DVectord &normals,
									 const DVectord &scale,
									 const std::string &name = "normal_field")
{
	using namespace libMesh;
	int mesh_dim = mesh.mesh_dimension();
	
	
	DenseVector<double> local_normal;
	DenseVector<Real> local_scale;
	
	std::vector<double> all_points, all_normals;
	
	std::vector<double> point(mesh_dim, 0.);
	for(auto n_it = mesh.active_nodes_begin(); n_it != mesh.active_nodes_end(); ++n_it) {
		Node &n = **n_it;
		
		std::vector<dof_id_type> node_dof_ids;
		
		for(int d = 0; d < mesh_dim; ++d) {
			auto dof_id = n.dof_number(0, 0, d);
			node_dof_ids.push_back(dof_id);
			
			point[d] = n(d);
		}
		
		get_vector(normals, node_dof_ids, local_normal);
		get_vector(scale,   node_dof_ids, local_scale);
		
		for(int d = 0; d < mesh_dim; ++d) {
			local_normal(d) *= local_scale(0);
		}
		
		if(local_normal.l2_norm() < 1e-16) continue;
		
		all_points.insert(all_points.end(),
						  point.begin(),
						  point.end());
		
		all_normals.insert(all_normals.end(),
						   local_normal.get_values().begin(),
						   local_normal.get_values().end());
	}
	
	if(all_points.empty()) {
		return;
	}
	
	quiver(mesh_dim,
		   all_points.size()/mesh_dim,
		   &all_points[0],
		   &all_normals[0],
		   name);
}


void run_biomechanics_example(libMesh::LibMeshInit &init)
{
	std::cout << "-----------------------------\n";
	std::cout << "run_biomechanics_example\n";
	//////////////////////////////////////////////////
	//////////////////////////////////////////////////
	
	// static const bool is_leaflet = true;
	// ContactSimParams params = leaflets_contact;
	
	//	static const bool is_leaflet = false;
	bool is_implant = false;
	// ContactSimParams params = contact_cylinder; static const int coords = 1;
	// ContactSimParams params = contact8; static const int coords = 1;
	// ContactSimParams params = contact_circles; static const int coords = 1;
	// ContactSimParams params = multi_contact_3D_2; static const int coords = 1;
	// ContactSimParams params = hip_femure_contact; static const int coords = 2;
	// ContactSimParams params = implant_contact; static const int coords = 1; is_implant = true;
	// ContactSimParams params = contact_cubes; static const int coords = 2;
	ContactSimParams params = hertz_contact; static const int coords = 1;
	//	ContactSimParams params = hertz_contact_coarse; static const int coords = 1;
	
	// predicate->add(101, 102);
	
	std::shared_ptr<cutlibpp::MasterAndSlave> predicate = nullptr;
	
	if(is_implant) {
		predicate = std::make_shared<cutlibpp::MasterAndSlave>();
		predicate->add(102, 101);
		predicate->add(103, 102);
		predicate->add(104, 103);
		predicate->add(105, 10);
	}
	
	auto mesh = make_shared<Mesh>(init.comm());
	// express::Communicator express_comm(init.comm().get());
	// mesh->read("/Users/patrick/Downloads/ASCII_bone/all_sidesets.e");
	
	// UGXMeshReader reader;
	// if(!reader.read("/Users/patrick/Downloads/AN_Keramik_Einlage_3971885250_3D01_96411_39-32.ugx", *mesh)) {
	// 	return;
	// }
	
	// MeshRefinement ref(*mesh);
	// ref.uniformly_coarsen();
	
	// plot_mesh(*mesh, "mesh");
	// return;
	
	mesh->read(params.mesh_path);
	// plot_mesh(*mesh, "mesh");
	
	const int dim = mesh->mesh_dimension();
	
	//////////////////////////////////////////////////
	//////////////////////////////////////////////////
	
	Chrono c;
	
	
	LibMeshFEContext<LinearImplicitSystem> context(mesh);
	auto space = vector_fe_space("stress_", LAGRANGE_VEC, FIRST, context);
	auto u = fe_function(space);
	
	//boundary conditions
	std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_1 = [dim, &params](const Point &, DenseVector<Real> &output) {
		output.resize(dim);
		output.zero();
		//offset to y coordinate of displacement
		output(coords) = params.dirichlet_value_1;
	};
	
	auto bc_1 = boundary_conditions(u == vec_coeff(boundary_cond_1),    { params.boundary_tag_1 });
	strong_enforce(bc_1);
	
	std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_2 = [dim, &params](const Point &, DenseVector<Real> &output) {
		output.resize(dim);
		output.zero();
		
		//offset to coordinate of displacement
		output(coords) = params.dirichlet_value_2;
	};
	
	auto bc_2 = boundary_conditions(u == vec_coeff(boundary_cond_2), { params.boundary_tag_2 });
	strong_enforce(bc_2);
	
	std::function<void (const Point &, DenseVector<Real> &output)> boundary_cond_3 = [dim, &params](const Point &, DenseVector<Real> &output) {
		output.resize(dim);
		output.zero();
		
		//offset to coordinate of displacement
		output(coords) = params.dirichlet_value_3;
		
	};
	
	if(params.boundary_tag_3 >= 0) {
		auto bc_3 = boundary_conditions(u == vec_coeff(boundary_cond_3), { params.boundary_tag_3 });
		strong_enforce(bc_3);
	}
	
	context.equation_systems.init();
	u.set_quad_rule(make_shared<libMesh::QGauss>(dim, SIXTH));
	
	
	DenseVector<Real> vec(dim);
	vec.zero();
	
	
	long n = u.dof_map().n_dofs();
	
	
	
	DSMatrixd coupling;
	DVectord gap;
	DVectord normals;
	DSMatrixd orhtogonal_trafos;
	
	
	
	std::vector<bool> is_contact_node;
	
	bool use_old_code = true;
	if(use_old_code) {
		MortarContactAssembler assembler(make_ref(space));
		assembler.set_strict_gap_policy(params.strict_search_radius);
		
		if(!assembler.assemble(coupling,
							   gap,
							   normals,
							   orhtogonal_trafos,
							   is_contact_node,
							   params.search_radius,
							   predicate
							   ))
		{
			//Just set some values that do not change the original system
			coupling = identity(n, n);
			orhtogonal_trafos = identity(n, n);
			gap = values(n, 10000);
			normals = zeros(n);
			is_contact_node.resize(n);
			std::fill(is_contact_node.begin(), is_contact_node.end(), false);
		}
	} else {
		//does not work with vector_fe_space
		//		DVectord is_contact_node_v, weighted_gap;
		//		std::shared_ptr<libMesh::DofMap> dof_map = utopia::make_ref(u.dof_map());
		//		DSMatrixd B;
		//
		//		express::Communicator express_comm(init.comm().get());
		//		utopia::assemble_contact(
		//								 express_comm,
		//								 mesh,
		//								 dof_map,
		//								 0,
		//								 B,
		//								 orhtogonal_trafos,
		//								 weighted_gap,
		//								 normals,
		//								 is_contact_node_v,
		//								 params.search_radius,
		//								 {{101, 102}},
		//								 true);
		//
		//		utopia::Range r = range(is_contact_node_v);
		//		is_contact_node.resize(r.extent(), false);
		//
		//		each_read(is_contact_node_v, [&is_contact_node, &r](const SizeType i, const double value) {
		//			is_contact_node[i - r.begin()] = value > 1e-16;
		//		});
		//
		//		DVectord d = sum(B, 1);
		//		auto s_B = local_size(B).get(0);
		//
		//		DSMatrixd D_inv = local_sparse(s_B, s_B, 1);
		//		{
		//			Write<DSMatrixd> v_d(D_inv);
		//			each_read(d, [&D_inv](const SizeType i, const double value) {
		//				if(value > 1e-15) {
		//					D_inv.set(i, i, 1./value);
		//				} else {
		//					D_inv.set(i, i, 1.);
		//				}
		//			});
		//		}
		//
		//		coupling = D_inv * B;
		//		coupling += local_identity(local_size(coupling));
		//		gap = D_inv * weighted_gap;
	}
	
	DSMatrixd K  = sparse(n, n, 20);
	DVectord rhs = zeros(n);
	
	double mu = 100.0, lambda = 200.0;
	LameeParameters lamee_params(mu, lambda);
	//	lamee_params.lambda_[1] = 5.;
	//	lamee_params.mu_[1] = 10.;
	
	c.start();
	auto ass = make_assembly([&]() -> void {
		assemble_linear_elasticity_system(lamee_params, u, *context.system.matrix, *context.system.rhs);
	});
	
	context.system.attach_assemble_object(ass);
	context.equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 1;
	context.equation_systems.solve();
	
	c.stop();
	std::cout << "Matrix and rhs assembly: ";
	c.describe(std::cout);
	std::cout << std::endl;
	
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
	newton.max_it(40);
	
	newton.set_box_constraints(make_upper_bound_constraints(make_ref(gap_c)));
	newton.solve(K_c, rhs_c, sol_c);
	
	// if(express_comm.isAlone()) plot_scaled_normal_field(*context.mesh, normals, gap_c);
	
	//Change back to original basis
	DVectord sol = coupling * (orhtogonal_trafos * sol_c);
	
	// convert(sol, *context.system.solution);
	
	plot_scaled_normal_field(*mesh, normals, gap, "biomech");
	
	DVectord stress;
	assemble_von_mises_stress(lamee_params, u, sol, stress);
	convert(stress, *context.system.solution);
	
//	{
//		DVectord m_contact_stress = (rhs - K * sol);
//		auto s_K = local_size(K);
//		DSMatrixd mass_matrix = local_sparse(s_K.get(0), s_K.get(1), 10);
//		assemble(u, u, integral(dot(u, u)), mass_matrix);
//		
//		stress = local_zeros(local_size(m_contact_stress));
//		solve(mass_matrix, m_contact_stress, stress);
//	}
	
	convert(stress, *context.system.solution);
	ExodusII_IO(*mesh).write_equation_systems ("elasticity_contact.e", context.equation_systems);
	
	int sys_num = u.dof_map().sys_number();
	int var_num = u.var_num();
	for(auto n_it = mesh->local_nodes_begin(); n_it != mesh->local_nodes_end(); ++n_it) {
		for(unsigned int c = 0; c != (*n_it)->n_comp(sys_num, var_num); ++c) {
			const int dof_id = (*n_it)->dof_number(sys_num, var_num, c);
			(**n_it)(c) += sol.get(dof_id);
		}
	}
	
	ExodusII_IO(*mesh).write_equation_systems ("elasticity_contact_deformed.e", context.equation_systems);
	
	
	
	// {
	// 	Read<DVectord> w_g(gap);
	// 	Read<DVectord> w_n(normals);
	
	// 	for(size_t i = 0; i < is_contact_node.size(); ++i) {
	// 		if(is_contact_node[i]) {
	// 			const double len = gap.get((i/dim)*dim);
	// 			context.system.solution->set(i, len * normals.get(i));
	// 		} else {
	// 			context.system.solution->set(i, 0);
	// 		}
	// 	}
	// }
	
	// ExodusII_IO(*mesh).write_equation_systems ("elasticity_contact_cn.e", context.equation_systems);
	
	
	
}

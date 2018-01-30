#include "utopia_FEEvalTest.hpp"

#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

#include "utopia_libmesh.hpp"
#include "libmesh/parallel_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/linear_implicit_system.h"

#include "utopia_LibMeshBackend.hpp"
#include "utopia_Equations.hpp"
#include "utopia_FEConstraints.hpp"
#include "utopia_FindSpace.hpp"
#include "utopia_IsForm.hpp"
#include "utopia_libmesh_NonLinearFEFunction.hpp"
#include "utopia_FEKernel.hpp"

#include "libmesh/exodusII_io.h"
#include <algorithm>

namespace utopia {

	template<class Derived, int Backend>
	static auto eval(const Expression<Derived> &expr, AssemblyContext<Backend> &ctx) 
	-> decltype( FEEval<Derived, Traits<Derived>, Backend, QUAD_DATA_NO>::apply(expr.derived(), ctx) )
	{
		return FEEval<Derived, Traits<Derived>, Backend, QUAD_DATA_NO>::apply(expr.derived(), ctx);
	}

	template<class Derived, int Backend>
	static auto quad_eval(const Expression<Derived> &expr, AssemblyContext<Backend> &ctx) 
	-> decltype( FEEval<Derived, Traits<Derived>, Backend, QUAD_DATA_YES>::apply(expr.derived(), ctx) )
	{
		return FEEval<Derived, Traits<Derived>, Backend, QUAD_DATA_YES>::apply(expr.derived(), ctx);
	}

	template<typename T>
	void disp(const libMesh::TensorValue<T> &t)
	{
		for(auto i = 0; i < LIBMESH_DIM; ++i) {
			for(auto j = 0; j < LIBMESH_DIM; ++j) {
				std::cout << t(i, j) << ", ";
			}

			std::cout << "\n";
		}
	}

	template<typename T>	
	void disp(const std::vector<std::vector<T>> &data) {
		for(std::size_t i = 0; i < data.size(); ++i) {
			for(std::size_t qp = 0; qp < data.size(); ++qp) {
				std::cout << "[" << i << ", " << qp << "]:\n";
				disp(data[i][qp]);
			}
		}
	}

	template<typename T>	
	void disp(const std::vector<T> &data) {
		for(std::size_t qp = 0; qp < data.size(); ++qp) {
			std::cout << "[" << qp << "]:\n";
			disp(data[qp]);
		}
	}

	libMesh::TensorValue<double> make_tensor_value(const LMDenseMatrix &in)
	{
		libMesh::TensorValue<double> ret;
		each_read(in, [&ret](const SizeType i, const SizeType j, const double val) {
			ret(i, j) = val;
		});

		return ret;
	}

	LMDenseMatrix make_wrapper(const int n, const libMesh::TensorValue<double> &in)
	{
		LMDenseMatrix ret = zeros(n, n);
		each_write(ret, [&in](const SizeType i, const SizeType j) ->  double {
			return in(i, j);
		});

		return ret;
	}


	LMDenseMatrix neohookean_first_piola(const double mu, const double lambda, const LMDenseMatrix &F_in)
	{
		libMesh::TensorValue<double> F = make_tensor_value(F_in);
		
		if(size(F_in).get(0) < 3) {
			F(2, 2) = 1;
		}

		libMesh::TensorValue<double> F_inv_t = F.inverse().transpose();
		const double J = F.det();

		libMesh::TensorValue<double>  P = mu * (F - F_inv_t) + (lambda * std::log(J)) * F_inv_t;
		return make_wrapper(size(F_in).get(0), P);
	}


	std::vector<LMDenseMatrix> neohookean_first_piola(const double mu, const double lambda, const std::vector<LMDenseMatrix> &F)
	{
		std::vector<LMDenseMatrix> ret = F;

		for(std::size_t i = 0; i < F.size(); ++i) {
			ret[i] = neohookean_first_piola(mu, lambda, F[i]);
		}

		return ret;
	}

	template<class Tensor>
	Tensor neohookean_linearized(const double mu, const double lambda, const Tensor &F)
	{
		Tensor F_inv_t = transpose(inv(F));
		const double J = det(F);
		return mu * (F - F_inv_t) + lambda * std::log(J) * F_inv_t;
	}

	void run_fe_eval_test(libMesh::LibMeshInit &init)
	{
		auto mesh = std::make_shared<libMesh::DistributedMesh>(init.comm());

		libMesh::MeshTools::Generation::build_square(*mesh,
			1, 1,
			0, 1,
			0, 1.,
			libMesh::QUAD4);

		auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);	
		auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("eval-test");

		const double mu = 1.;
		const double lambda = 1.;

		auto elem_order = libMesh::FIRST;

		////////////////////////////////////////////

		auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
		auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
		auto V = Vx * Vy;

		auto u = trial(V);
		auto v = test(V);

		Vx.initialize();
		DVectord sol = ghosted(Vx.dof_map().n_local_dofs(), Vx.dof_map().n_dofs(), Vx.dof_map().get_send_list());
		sol = local_values(local_size(sol).get(0), 1.);
		
		{
			Write<DVectord> w_sol(sol);
			sol.set(range(sol).begin(), .0);
		}

		//
		AssemblyContext<LIBMESH_TAG> ctx;
		ctx.set_current_element(0);
		ctx.set_has_assembled(false);
		ctx.init_bilinear( inner(u, v) * dX + inner(grad(u), grad(v)) * dX );

		//symbolic expressions
		auto uk 	 = interpolate(sol, u);
		auto g_uk    = grad(uk);
		auto F 		 = identity() + g_uk;
		auto F_inv   = inv(F);
		auto F_inv_t = transpose(F_inv);
		auto J 		 = det(F);

		auto P = mu * (F - F_inv_t) + (lambda * logn(J)) * F_inv_t;

		auto stress_lin = mu * grad(u) 
		-(lambda * logn(J) - mu) * F_inv_t * transpose(grad(u)) * F_inv_t 
		+ inner(lambda * F_inv_t, grad(u)) * F_inv_t;

		//evaluate 
		auto eval_uk 	  = eval(uk, ctx);
		auto eval_g_uk	  = eval(g_uk, ctx);
		auto eval_F 	  = eval(F, ctx);
		auto eval_F_inv   = eval(F_inv, ctx);
		auto eval_F_inv_t = eval(F_inv_t, ctx);
		auto eval_J       = eval(J, ctx);

		
		auto eval_P = quad_eval(P, ctx);
		auto eval_P_expected = neohookean_first_piola(mu, lambda, eval_F);

		for(std::size_t i = 0; i < eval_P.size(); ++i) {
			libMesh::DenseMatrix<double> mat = eval_P[i].implementation();
			mat.add(-1., eval_P_expected[i].implementation());
			assert(approxeq(mat.linfty_norm(), 0.));
		}

		auto eval_stress = quad_eval(stress_lin, ctx);

		disp(eval_stress);

	}
}
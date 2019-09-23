#include "utopia_libmesh.hpp"
#include "utopia_FEEvalTest.hpp"

#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

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

    // LMDenseMatrix make_tensor_value(const LMDenseMatrix &in)
    // {
    //     LMDenseMatrix ret;
    //     each_read(in, [&ret](const SizeType i, const SizeType j, const double val) {
    //         ret(i, j) = val;
    //     });

    //     return ret;
    // }

    // LMDenseMatrix make_wrapper(const int n, const LMDenseMatrix &in)
    // {
    //     LMDenseMatrix ret = zeros(n, n);
    //     each_write(ret, [&in](const SizeType i, const SizeType j) ->  double {
    //         return in(i, j);
    //     });

    //     return ret;
    // }

    LMDenseMatrix neohookean_first_piola(const double mu, const double lambda, const LMDenseMatrix &F)
    {
        LMDenseMatrix F_inv_t = transpose(inv(F));
        return mu * (F - F_inv_t) + (lambda * std::log(det(F))) * F_inv_t;
    }

    // LMDenseMatrix neohookean_first_piola(const double mu, const double lambda, const LMDenseMatrix &F)
    // {
    //     // LMDenseMatrix F = F;
    //     // if(size(F).get(0) < 3) {
    //     //     F.set(2, 2, 1);
    //     // }

    //     return neohookean_first_piola(mu, lambda, F);
    // }

    std::vector<LMDenseMatrix> neohookean_first_piola(const double mu, const double lambda, const std::vector<LMDenseMatrix> &F)
    {
        std::vector<LMDenseMatrix> ret = F;

        for(std::size_t i = 0; i < F.size(); ++i) {
            ret[i] = neohookean_first_piola(mu, lambda, F[i]);
        }

        return ret;
    }

    template<class Tensor>
    Tensor neohookean_linearized(const double mu, const double lambda, const Tensor &H, const Tensor &F)
    {
        Tensor F_inv_t = transpose(inv(F));
        const double J = det(F);
        const double alpha = (1.0 * lambda * std::log(J) - 1.0 * mu);

        return mu * H - alpha * F_inv_t * transpose(H) * F_inv_t + lambda * inner(F_inv_t, H) * F_inv_t;
    }

    std::vector<std::vector<LMDenseMatrix>> neohookean_linearized(const double mu, const double lambda, const std::vector<std::vector<LMDenseMatrix>> &H, const std::vector<LMDenseMatrix> &F)
    {
        std::vector<std::vector<LMDenseMatrix>> ret;

        ret.resize(H.size());

        const auto dim = size(F[0]).get(0);

        for(std::size_t i = 0; i < H.size(); ++i) {
            ret[i].resize(F.size());

            assert(H[i].size() == F.size());

            for(std::size_t qp = 0; qp < F.size(); ++qp) {


                // if(dim < 3) {
                //     F_qp(2, 2) = 1;
                // }

                auto ret_iqp = neohookean_linearized(mu, lambda, H[i][qp], F[qp]);

                // if(dim < 3) {
                //     ret_iqp.set(2, 2, 0;
                //     ret_iqp(2, 0) = 0;
                //     ret_iqp(2, 1) = 0;
                //     ret_iqp(0, 2) = 0;
                //     ret_iqp(1, 2) = 0;
                // }

                ret[i][qp] = ret_iqp;
            }
        }

        return ret;
    }

    template<class T>
    void check_equal(const std::vector<T> &left, const std::vector<T> &right)
    {
        for(std::size_t i = 0; i < left.size(); ++i) {
            assert(approxeq(left[i], right[i], 0.));
        }
    }

    void check_equal(
        const std::vector<std::vector<LMDenseMatrix>> &left,
        const std::vector<std::vector<LMDenseMatrix>> &right)
    {
        assert(left.size() == right.size());

        for(std::size_t i = 0; i < left.size(); ++i) {
            assert(left[i].size() == right[i].size());

            for(std::size_t qp = 0; qp < left[i].size(); ++qp) {
                auto l = left[i][qp];
                auto r = right[i][qp];

                const bool ok = approxeq(l, r, 0.);

                if(!ok) {
                    disp("----------------");
                    disp(l);
                    disp("----------------");
                    disp(r);
                    disp("----------------");

                }
                assert(ok);
            }
        }
    }

    void FEEvalTest::run(Input &in)
    {
        auto mesh = std::make_shared<libMesh::DistributedMesh>(comm());

        const int n = 1;
        libMesh::MeshTools::Generation::build_square(*mesh,
            n, n,
            0, 1,
            0, 1.,
            // libMesh::QUAD4
            libMesh::TRI3
            );

        auto dim = mesh->mesh_dimension();

        auto equation_systems = std::make_shared<libMesh::EquationSystems>(*mesh);
        auto &sys = equation_systems->add_system<libMesh::LinearImplicitSystem>("eval-test");

        const double mu = 0.1;
        const double lambda = 0.1;

        auto elem_order = libMesh::FIRST;

        ////////////////////////////////////////////

        auto Vx = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_x");
        auto Vy = LibMeshFunctionSpace(equation_systems, libMesh::LAGRANGE, elem_order, "disp_y");
        auto V = Vx * Vy;

        auto u = trial(V);
        auto v = test(V);

        Vx.initialize();

        UIndexSet ghost_nodes;
        convert(Vx.dof_map().get_send_list(), ghost_nodes);
        UVector sol = ghosted(Vx.dof_map().n_local_dofs(), Vx.dof_map().n_dofs(), ghost_nodes);
        sol = local_values(local_size(sol).get(0), 0.);

        {
            Write<UVector> w_sol(sol);
            auto r = range(sol);
            // for(auto i = r.begin(); i != r.end(); ++i) {
            // 	sol.set(i, i + (i%3)*1.5);
            // }

            sol.set(2, .5);
        }

        disp(sol);

        libMesh::ExodusII_IO io(*mesh);
        convert(sol, *sys.solution);
        sys.solution->close();
        io.write_equation_systems("m.e", *equation_systems);

        AssemblyContext<LIBMESH_TAG> ctx;
        for(auto e_it = mesh->active_local_elements_begin(); e_it != mesh->active_local_elements_end(); ++e_it) {
            //NEOHOOKEAN---------------------------------------------------
            //symbolic expressions
            auto uk 	 = interpolate(sol, u);
            auto g_uk    = grad(uk);
            auto F 		 = identity() + g_uk;
            auto F_t 	 = transpose(F);
            auto F_inv   = inv(F);
            auto F_inv_t = transpose(F_inv);
            auto J 		 = det(F);

            auto P = mu * (F - F_inv_t) + (lambda * logn(J)) * F_inv_t;

            auto mixedUp = mu * grad(u) - F; //(lambda * logn(J) - mu) * F_inv_t * transpose(grad(u)) * F_inv_t;

            auto stress_lin = mu * grad(u)
            -(lambda * logn(J) - mu) * F_inv_t * transpose(grad(u)) * F_inv_t
            + inner(lambda * F_inv_t, grad(u)) * F_inv_t;

            ctx.set_current_element((*e_it)->id());
            ctx.set_has_assembled(false);
            ctx.init( inner(stress_lin, grad(v)) * dX == inner(P, grad(v)) );

            //evaluate
            auto eval_grad    = eval(grad(u), ctx);
            auto eval_uk 	  = eval(uk, ctx);
            auto eval_g_uk	  = eval(g_uk, ctx);
            auto eval_F 	  = eval(F, ctx);
            auto eval_F_inv   = eval(F_inv, ctx);
            auto eval_F_inv_t = eval(F_inv_t, ctx);
            auto eval_J       = eval(J, ctx);
            auto eval_log_J   = eval(logn(J), ctx);

            auto eval_mixedUp = quad_eval(mixedUp, ctx);

            disp("-----------------------------------");
            disp("-----------------------------------");
            // disp(eval_uk);
            // disp("-----------------------------------");
            // disp(eval_grad);
            // disp("-----------------------------------");
            // disp(eval_g_uk);
            // disp("-----------------------------------");
            // disp(eval_F);
            // disp("-----------------------------------");

            // disp(eval_F_inv);
            // disp("-----------------------------------");

            auto eval_P = quad_eval(P, ctx);
            auto eval_P_expected = neohookean_first_piola(mu, lambda, eval_F);

            check_equal(eval_P, eval_P_expected);

            // auto eval_stress = quad_eval(stress_lin, ctx);
            // auto eval_stress_expected = neohookean_linearized(mu, lambda, eval_grad, eval_F);
            // check_equal(eval_stress, eval_stress_expected);


            // auto g = eval(inner(P, grad(v)) * dX, ctx);
            // // disp(g);
            // disp("-----------------------------------");
            // auto H = eval(inner(stress_lin, grad(v)) * dX, ctx);
            // // disp(H);


            // //////////////////////////////////////////////////////////////

            // auto C = F_t * F;
            // auto E = 0.5 * (C - identity());
            // auto S = 2.0 * mu * E + lambda * (trace(E) * identity());
            // auto C_lin = 0.5 * (F_t * grad(u) + transpose(grad(u)) * F);

            // auto stress_lin_2 = F * ( (2.0 * mu) * C_lin + lambda * (trace(C_lin) * identity()) ) + grad(u) * S;

            // auto eval_C = eval(C, ctx);
            // auto eval_E = eval(E, ctx);
            // auto eval_trace_ExI = eval(trace(E) * identity(), ctx);
            // auto eval_S = eval(S, ctx);
            // auto eval_C_lin = eval(C_lin, ctx);
            // auto eval_stress_lin_2 = eval(stress_lin_2, ctx);


            // ///////////////////////////////////////////////////////////////////

            // auto S_bar = mu * identity(dim, dim);
            // auto eval_S_bar = quad_eval(S_bar * F, ctx);
            // // disp(eval_S_bar);

            // // auto sum_SF = F - S_bar;
            // auto sum_SF = S_bar + F;
            // // auto sum_SF = mu * identity() - F;
            // auto eval_sum_SF = eval(inner(sum_SF, grad(v)), ctx);
            // // disp(eval_sum_SF);

            // auto S_iso = S_bar + (inner((-1.0 / 2) * S_bar, C) * inv(C));
            // auto eval_S_iso = quad_eval(S_iso, ctx);
            // // MostDescriptive<decltype(S_bar), decltype(F)>::Type desc;
            // // std::cout <<

            // auto dot_grads = inner(grad(uk), grad(uk));
            // auto en = quad_eval(dot_grads, ctx);
            // assert(!en.empty());

            // auto dot_grads_dx = dot_grads * dX;
            // auto endx = eval(dot_grads_dx, ctx);
            // // disp(endx);

            // auto denom      = inner(grad(uk), grad(uk));
            // auto div_inner  = inner(grad(uk)/denom, grad(v));

            // auto e_div_inner = quad_eval(div_inner, ctx);
            // disp(e_div_inner[0]);

            // auto div_inner_dx = div_inner * dX;

            // auto e_div_inner_dx = eval(div_inner_dx, ctx);
            // disp(e_div_inner_dx);
        }
    }
}
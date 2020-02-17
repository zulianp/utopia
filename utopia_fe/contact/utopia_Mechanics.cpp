#include "utopia_Mechanics.hpp"
#include "utopia_Contact.hpp"
#include "libmesh/dof_map.h"

//REMOVE ME
#include "utopia_LibMeshBackend.hpp"

#include "utopia_FormEvaluator.hpp"
#include "utopia_fe_core.hpp"
#include "utopia.hpp"
// //#include "utopia_fe_homemade.hpp"
#include "utopia_FEIsSubTree.hpp"
#include "utopia_MixedFunctionSpace.hpp"

#include "utopia_libmesh_NonLinearFEFunction.hpp"

#include "utopia_libmesh.hpp"
#include <cmath>

namespace utopia {

    void MechanicsContext::init_mass_matrix(const ProductFunctionSpace<LibMeshFunctionSpace> &V)
    {
        auto u = trial(V);
        auto v = test(V);

        assemble(inner(u, v) * dX, non_lumped_mass_matrix);
        UVector lumped_mass_vector = sum(non_lumped_mass_matrix, 1);
        mass_matrix = diag(lumped_mass_vector);
        inverse_mass_vector = 1./lumped_mass_vector;
    }

    void MechanicsState::init(const Size &local_size, const Size &global_size)
    {
        displacement = local_zeros(local_size);
        displacement_increment = local_zeros(local_size);

        velocity = local_zeros(local_size);

        internal_force = local_zeros(local_size);
        external_force = local_zeros(local_size);

        stress = local_zeros(local_size);
        t = 0.;
    }


    MechWithContactIntegrationScheme::MechWithContactIntegrationScheme(
        const unsigned int dim,
        libMesh::DofMap &dof_map)
    : dim(dim), dof_map(dof_map)
    {
        linear_solver = std::make_shared<Factorization<USparseMatrix, UVector>>();
        // auto temp =  std::make_shared<BiCGStab<USparseMatrix, UVector>>();
        // auto temp =  std::make_shared<ConjugateGradient<USparseMatrix, UVector, HOMEMADE>>();
        // auto temp =  std::make_shared<GaussSeidel<USparseMatrix, UVector>>();

        // temp->atol(1e-10);
        // temp->stol(1e-10);
        // temp->rtol(1e-10);
        // // temp->verbose(true);
        // temp->max_it(std::max(long(2), long(dof_map.n_dofs()/2)));

        // linear_solver =  temp;
    }

    bool MechWithContactIntegrationScheme::solve(
        const USparseMatrix &K,
        const UVector &inverse_mass_vector,
        const UVector &rhs,
        const UVector &gap,
        const Friction &friction,
        UVector &sol)
    {
        //FIXME
        UVector non_const_gap = gap;

        if(friction.friction_coefficient == 0.) {
            if(mpi_world_size() == 1) {
                SemismoothNewton<USparseMatrix, UVector> solver(linear_solver);
                solver.max_it(100);
                solver.verbose(true);
                solver.set_box_constraints(make_upper_bound_constraints(make_ref(non_const_gap)));
                return solver.solve(K, rhs, sol);

            } else {
                // SemismoothNewton<USparseMatrix, UVector, PETSC_EXPERIMENTAL> solver(linear_solver);
                SemismoothNewton<USparseMatrix, UVector> solver(linear_solver);
                solver.max_it(100);

                // ProjectedGaussSeidel<USparseMatrix, UVector, -1> solver;
                // solver.max_it(size(rhs).get(0) * 40);

                solver.verbose(true);
                solver.set_box_constraints(make_upper_bound_constraints(make_ref(non_const_gap)));
                return solver.solve(K, rhs, sol);
            }

        } else {
            typedef std::function<void(const USparseMatrix &, const UVector &, const UVector &, UVector &, UVector &)> F;

            UVector lambda, d;
            F f = [&](const USparseMatrix &H, const UVector &g, const UVector &x, UVector &active, UVector &value)
            {
                lambda = (g - H * x);
                lambda = e_mul(inverse_mass_vector, lambda);
                apply_zero_boundary_conditions(dof_map, lambda);

                d = lambda + (x - non_const_gap);

                Read<UVector> r_d(d);
                Read<UVector> r_l(lambda);
                Read<UVector> r_u(non_const_gap);

                Write<UVector> w_d(active);
                Write<UVector> w_v(value);

                auto rr = range(x);
                for (SizeType i = rr.begin(); i != rr.end(); i += dim) {
                    if (d.get(i) >= -1e-16) {
                        active.set(i, 1.0);
                        value.set(i, non_const_gap.get(i));

                        double n_s = lambda.get(i);
                        double t_s = 0.;

                        for(SizeType d = 1; d < dim; ++d) {
                            const double t_d = lambda.get(i + d);
                            t_s = t_d * t_d;
                        }

                        t_s = std::sqrt(t_s);

                        if(t_s < friction.friction_coefficient * n_s) {
                            for(SizeType d = 1; d < dim; ++d) {
                                active.set(i + d, 1.0);
                                value.set(i + d,  0.0);
                            }

                        } else {

                            for(SizeType d = 1; d < dim; ++d) {
                                active.set(i + d, 0.0);
                                value.set(i + d,  0.0);
                            }
                        }

                    } else {
                        for(SizeType d = 0; d < dim; ++d) {
                            active.set(i + d, 0.0);
                            value.set(i + d,  0.);
                        }
                    }
                }
            };



            GenericSemismoothNewton<USparseMatrix, UVector, F> solver(f, linear_solver);
            solver.verbose(true);
            solver.max_it(20);
            return solver.solve(K, rhs, sol);
        }
    }

    void ImplicitEuler::apply(
        const double dt,
        const MechanicsContext &mech_ctx,
        const MechanicsState &old,
        MechanicsState &current)
    {
        const USparseMatrix &K = mech_ctx.stiffness_matrix;
        UVector rhs = current.external_force - old.internal_force;

        linear_solver->solve(K, rhs, current.displacement_increment);
        current.displacement = old.displacement + current.displacement_increment;
        current.internal_force = K * current.displacement;
        current.t = old.t + dt;

        //FIXME find other way
        apply_zero_boundary_conditions(dof_map, current.internal_force);
    }

    void ImplicitEuler::apply(
        const double dt,
        const MechanicsContext &mech_ctx,
        const Contact  &contact,
        const Friction &friction,
        const MechanicsState &old,
        MechanicsState &current)
    {
        auto s = local_size(old.internal_force);

        // const USparseMatrix &T = contact.complete_transformation;
        const USparseMatrix &K = mech_ctx.stiffness_matrix;
        const UVector rhs = current.external_force;// - old.internal_force;

        UVector sol_c = local_zeros(s);
        UVector rhs_c;// = transpose(T) * rhs;
        USparseMatrix K_c;//  = transpose(T) * K * T;

        contact.couple(rhs, rhs_c);
        contact.couple(K, K_c);

        bool solved = solve(K_c, mech_ctx.inverse_mass_vector, rhs_c, contact.gap(), friction, sol_c);

        if(!solved) {
            std::cerr << "[Error] unable to solve non-linear system" << std::endl;
        }

        assert(solved);

        // current.displacement_increment = T * sol_c;
        contact.uncouple(sol_c, current.displacement_increment);
        current.displacement = old.displacement + current.displacement_increment;
        current.internal_force = K * current.displacement_increment;
        current.t = old.t + dt;

        // current.stress = e_mul(mech_ctx.inverse_mass_vector, T * (rhs_c - K_c * sol_c));
        // current.stress = T * (rhs_c - K_c * sol_c);

        // current.stress = e_mul(contact.inv_mass_vector, (current.external_force - current.internal_force));
        // current.stress = e_mul(contact.inv_mass_vector(), (current.external_force - current.internal_force));
        contact.remove_mass((current.external_force - current.internal_force), current.stress);
        //FIXME find other way
        apply_zero_boundary_conditions(dof_map, current.internal_force);
        apply_zero_boundary_conditions(dof_map, current.stress);
    }
}

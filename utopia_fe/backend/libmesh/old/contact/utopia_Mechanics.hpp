#ifndef UTOPIA_MECHANICS_HPP
#define UTOPIA_MECHANICS_HPP

#include "utopia.hpp"
#include "utopia_FEForwardDeclarations.hpp"
#include "utopia_fe_base.hpp"
#include "utopia_libmesh_FEForwardDeclarations.hpp"
#include "utopia_libmesh_Types.hpp"

#include <memory>

namespace libMesh {
    class DofMap;
}

namespace utopia {
    class Contact;

    class MechanicsState {
    public:
        void init(const Size &local_size, const Size &global_size);

        UVector displacement;
        UVector displacement_increment;

        UVector velocity;

        UVector internal_force;
        UVector external_force;

        UVector stress;

        double t;
    };

    class MechanicsContext {
    public:
        // stiffness matrix
        USparseMatrix stiffness_matrix;

        USparseMatrix non_lumped_mass_matrix;
        // lumped mass matrix
        USparseMatrix mass_matrix;
        UVector inverse_mass_vector;

        UVector dirichlet_selector;

        void init_mass_matrix(const ProductFunctionSpace<LibMeshFunctionSpace> &V);
    };

    class Friction {
    public:
        Friction() : friction_coefficient(0.) {}

        // FIXME
        double friction_coefficient;
    };

    class MechIntegrationScheme {
    public:
        virtual ~MechIntegrationScheme() {}

        virtual void apply(const double dt,
                           const MechanicsContext &mech_ctx,
                           const MechanicsState &old,
                           MechanicsState &current) = 0;
    };

    class MechWithContactIntegrationScheme {
    public:
        virtual ~MechWithContactIntegrationScheme() {}

        MechWithContactIntegrationScheme(const unsigned int dim, libMesh::DofMap &dof_map);

        void set_linear_solver(const std::shared_ptr<LinearSolver<USparseMatrix, UVector> > &solver) {
            this->linear_solver = solver;
        }

        virtual void apply(const double dt,
                           const MechanicsContext &mech_ctx,
                           const Contact &contact,
                           const Friction &friction,
                           const MechanicsState &old,
                           MechanicsState &current) = 0;

        bool solve(const USparseMatrix &K,
                   const UVector &inverse_mass_vector,
                   const UVector &rhs,
                   const UVector &gap,
                   const Friction &friction,
                   UVector &sol);

        // FIXME I do not like it
        unsigned int dim;
        libMesh::DofMap &dof_map;
        std::shared_ptr<LinearSolver<USparseMatrix, UVector> > linear_solver;
    };

    class ImplicitEuler : public MechIntegrationScheme, public MechWithContactIntegrationScheme {
    public:
        ImplicitEuler(const unsigned int dim, libMesh::DofMap &dof_map)
            : MechWithContactIntegrationScheme(dim, dof_map) {}

        void apply(const double dt,
                   const MechanicsContext &mech_ctx,
                   const MechanicsState &old,
                   MechanicsState &current) override;

        void apply(const double dt,
                   const MechanicsContext &mech_ctx,
                   const Contact &contact,
                   const Friction &friction,
                   const MechanicsState &old,
                   MechanicsState &current) override;
    };

    class ExternalForce {
    public:
        virtual ~ExternalForce() {}
        virtual void eval(const double t, UVector &result) = 0;
    };

    class ConstantExternalForce : public ExternalForce {
    public:
        inline void eval(const double, UVector &result) override { result = value; }

        template <class LinearForm>
        void init(const LinearForm &linear_form) {
            assemble(linear_form, value);
        }

        UVector value;
    };

    // class MechanicsSimulation {
    // public:

    // };
}  // namespace utopia

#endif  // UTOPIA_MECHANICS_HPP

#ifndef UTOPIA_WEAR_HPP
#define UTOPIA_WEAR_HPP

#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"

#include "utopia_FEForwardDeclarations.hpp"

#include "utopia.hpp"
#include "utopia_ContactAssembler.hpp"
#include "utopia_ContactSolver.hpp"
#include "utopia_Mechanics.hpp"
#include "utopia_fe_EDSL.hpp"
#include "utopia_libmesh_old.hpp"

#include <cmath>
#include <vector>

namespace utopia {

    void apply_displacement(const UVector &displacement_increment,
                            const libMesh::DofMap &dof_map,
                            libMesh::MeshBase &mesh);

    class Wear {
    public:
        using ContactSolverT = utopia::ContactSolver<USparseMatrix, UVector>;
        using ContactT = ContactSolverT::ContactT;

        Wear();

        void set_params(const double wear_coefficient, const double extrapolation_factor) {
            this->wear_coefficient = wear_coefficient;
            this->extrapolation_factor = extrapolation_factor;
        }

        inline void update(const double dt, const UVector &sliding_distance, const UVector &normal_stress) {
            if (empty(wear)) {
                wear.zeros(layout(normal_stress));
            }

            wear += (dt * wear_coefficient) * abs(e_mul(sliding_distance, normal_stress));
        }

        void compute_displacement(ProductFunctionSpace<LibMeshFunctionSpace> &V,
                                  const std::vector<int> &boundary_tags,
                                  UVector &wear_induced_displacement);

        void mesh_displacement(ProductFunctionSpace<LibMeshFunctionSpace> &V,
                               const std::vector<int> &boundary_tags,
                               UVector &disp);

        void init_aux_system(libMesh::EquationSystems &es, libMesh::Order order = libMesh::FIRST);

        void update_aux_system(const int main_system_number,
                               const MechanicsState &state,
                               const ContactT &contact,
                               const double dt,
                               libMesh::EquationSystems &es);

        void print(std::ostream &os) {
            for (auto w : total_wear) {
                os << w << " ";
            }

            os << "\n";
        }

    private:
        UVector wear;
        double wear_coefficient;
        double extrapolation_factor;

        std::vector<int> var_num_aux;
        unsigned int param_sys_number;

        std::vector<double> total_wear;

        // buffers
        UVector wear_induced_displacement;
        UVector is_normal_component;
        UVector normals;
        USparseMatrix trafo;
    };
}  // namespace utopia

#endif  // UTOPIA_WEAR_HPP

#ifndef UTOPIA_FRACTURE_FLOW_UTILS_HPP
#define UTOPIA_FRACTURE_FLOW_UTILS_HPP

#include "utopia.hpp"
#include "utopia_SemiGeometricMultigrid.hpp"

#include "libmesh/mesh.h"
#include "libmesh/nemesis_io.h"

#include <memory>

namespace utopia {

    std::shared_ptr<SemiGeometricMultigrid> make_mg_solver(const LibMeshFunctionSpace &space, const int n_levels);

    void write_solution(const std::string &name,
                        UVector &sol,
                        const int time_step,
                        const double t,
                        LibMeshFunctionSpace &space,
                        libMesh::Nemesis_IO &io);

    void transform_values(const LibMeshFunctionSpace &from,
                          const UVector &from_vector,
                          const LibMeshFunctionSpace &to,
                          UVector &to_vector,
                          std::function<double(const double &)> fun);

    void copy_values(const LibMeshFunctionSpace &from,
                     const UVector &from_vector,
                     const LibMeshFunctionSpace &to,
                     UVector &to_vector);

    void remove_constrained_dofs(libMesh::DofMap &dof_map, USparseMatrix &mat);

}  // namespace utopia

#endif  // UTOPIA_FRACTURE_FLOW_UTILS_HPP

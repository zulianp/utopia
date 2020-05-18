#include "utopia_FractureFlowUtils.hpp"
#include "utopia_InterpolationLocalAssembler.hpp"
#include "utopia_TransferAssembler.hpp"
#include "utopia_assemble_volume_transfer.hpp"

#include "moonolith_communicator.hpp"

#include "libmesh/elem.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"

#include <vector>

namespace utopia {

    std::shared_ptr<SemiGeometricMultigrid> make_mg_solver(const LibMeshFunctionSpace &space, const int n_levels) {
        auto linear_solver = std::make_shared<Factorization<USparseMatrix, UVector>>();
        auto smoother = std::make_shared<GaussSeidel<USparseMatrix, UVector>>();
        auto mg = std::make_shared<SemiGeometricMultigrid>(smoother, linear_solver);

        mg->algebraic().rtol(1e-9);
        mg->algebraic().atol(1e-14);
        // mg->verbose(true);
        mg->init(space, n_levels);
        return mg;
    }

    void write_solution(const std::string &name,
                        UVector &sol,
                        const int time_step,
                        const double t,
                        LibMeshFunctionSpace &space,
                        libMesh::Nemesis_IO &io) {
        utopia::convert(sol, *space.equation_system().solution);
        space.equation_system().solution->close();
        io.write_timestep(name, space.equation_systems(), time_step, t);
    }

    void transform_values(const LibMeshFunctionSpace &from,
                          const UVector &from_vector,
                          const LibMeshFunctionSpace &to,
                          UVector &to_vector,
                          std::function<double(const double &)> fun) {
        using libMesh::dof_id_type;

        if (empty(to_vector)) {
            to_vector = local_zeros(to.dof_map().n_local_dofs());
        }

        const auto &m = from.mesh();
        assert(&m == &to.mesh());

        std::vector<dof_id_type> from_dofs, to_dofs;

        Read<UVector> r_(from_vector);
        Write<UVector> w_(to_vector);

        auto r = range(to_vector);

        for (auto e_it = elements_begin(m); e_it != elements_end(m); ++e_it) {
            const auto &e = *e_it;

            from.dof_map().dof_indices(e, from_dofs, from.subspace_id());
            to.dof_map().dof_indices(e, to_dofs, to.subspace_id());

            assert(from_dofs.size() == to_dofs.size());

            auto n = from_dofs.size();

            for (std::size_t i = 0; i < n; ++i) {
                if (r.inside(to_dofs[i])) {
                    to_vector.set(to_dofs[i], fun(from_vector.get(from_dofs[i])));
                }
            }
        }
    }

    void copy_values(const LibMeshFunctionSpace &from,
                     const UVector &from_vector,
                     const LibMeshFunctionSpace &to,
                     UVector &to_vector) {
        transform_values(from, from_vector, to, to_vector, [](const double &value) -> double { return value; });
        // using libMesh::dof_id_type;

        // if(empty(to_vector)) {
        // 	to_vector = local_zeros(to.dof_map().n_local_dofs());
        // }

        // const auto &m = from.mesh();
        // assert(&m == &to.mesh());

        // std::vector<dof_id_type> from_dofs, to_dofs;

        // Read<UVector> r_(from_vector);
        // Write<UVector> w_(to_vector);

        // for(auto e_it = elements_begin(m); e_it != elements_end(m); ++e_it) {
        // 	const auto &e = *e_it;

        // 	from.dof_map().dof_indices(e, from_dofs, from.subspace_id());
        // 	to.dof_map().dof_indices(e, to_dofs, to.subspace_id());

        // 	assert(from_dofs.size() == to_dofs.size());

        // 	auto n = from_dofs.size();

        // 	for(std::size_t i = 0; i < n; ++i) {
        // 		to_vector.set(to_dofs[i], from_vector.get(from_dofs[i]));
        // 	}
        // }
    }

    void remove_constrained_dofs(libMesh::DofMap &dof_map, USparseMatrix &mat) {
        if (utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions begin: " << std::endl;
        }

        Chrono c;
        c.start();

        assert(!empty(mat));

        using SizeType = Traits<UVector>::SizeType;

        const bool has_constaints = dof_map.constraint_rows_begin() != dof_map.constraint_rows_end();

        Size ls = local_size(mat);
        Size s = size(mat);

        std::vector<SizeType> index;

        Range rr = row_range(mat);

        if (has_constaints) {
            for (SizeType i = rr.begin(); i < rr.end(); ++i) {
                if (dof_map.is_constrained_dof(i)) {
                    index.push_back(i);
                }
            }
        }

        set_zero_rows(mat, index, 0.);

        c.stop();

        if (utopia::Utopia::instance().verbose()) {
            std::cout << "apply_boundary_conditions end: " << c << std::endl;
        }
    }

}  // namespace utopia

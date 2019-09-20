#include "utopia_SemiGeometricMultigrid.hpp"

#include "utopia_libmesh.hpp"
#include "utopia_Socket.hpp"
#include "utopia_assemble_volume_transfer.hpp"
#include "utopia_Contact.hpp"
#include "utopia_BoundingBoxCoarsener.hpp"
#include "utopia_make_unique.hpp"

#include "moonolith_communicator.hpp"



#include "libmesh/mesh_tools.h"
#include "libmesh/explicit_system.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/libmesh_version.h"

#include <cmath>

namespace utopia {
    inline double volume(LibMeshFunctionSpace &V, int block = -1)
    {
        auto &dof_map = V.dof_map();
        auto u = trial(V);

        Traits<UVector>::IndexArray ghost_nodes;
        convert(dof_map.get_send_list(), ghost_nodes);
        UVector x = ghosted(dof_map.n_local_dofs(), dof_map.n_dofs(), ghost_nodes);
        x.set(1.);

        double volume = -1.;

        utopia::assemble(
            integral(interpolate(x, u), block),
            volume
        );

        return volume;
    }

    static void make_d(const USparseMatrix &mat, UVector &res)
    {
        res = sum(mat, 1);


        // std::cout << "make_d sum(mat): " << double(sum(res)) << std::endl;

        ReadAndWrite<UVector> rw_(res);
        auto r = range(res);
        for(auto k = r.begin(); k != r.end(); ++k) {
            if(approxeq(res.get(k), 0.0, 1e-14)) {
                res.set(k, 1.);
            }
        }
    }

    SemiGeometricMultigrid::SemiGeometricMultigrid(
        const std::shared_ptr<Smoother<USparseMatrix, UVector> > &smoother,
        const std::shared_ptr<LinearSolver<USparseMatrix, UVector> > &linear_solver)
    : mg(smoother, linear_solver),
      is_block_solver_(false),
      separate_subdomains_(false),
      use_interpolation_(false),
      use_coarse_interpolators_(true)
    { }

    void SemiGeometricMultigrid::generate_coarse_meshes(
        const libMesh::MeshBase &fine_mesh,
        const std::size_t n_levels,
        const int order_fine_level)
    {
        assert(order_fine_level == 1 || order_fine_level == 2);

        const std::size_t n_coarse_spaces = n_levels - 1;
        bool is_second_order = 2 == order_fine_level;

        meshes.resize(n_coarse_spaces);
        equation_systems.resize(n_coarse_spaces);
        use_interpolation_at_level_.resize(n_coarse_spaces);
        std::fill(std::begin(use_interpolation_at_level_), std::end(use_interpolation_at_level_), use_interpolation_);

        if(order_fine_level > 1) {
            meshes[n_coarse_spaces - 1] = fine_mesh.clone();
            use_interpolation_at_level_[0] = true;
            if(n_levels == 2) {
                return;
            }
        }

        const std::size_t n_p1_levels = n_coarse_spaces - is_second_order;
        // meshes[0] = generate_box_mesh(fine_mesh, n_levels - is_second_order);
        meshes[0] = generate_box_mesh(fine_mesh, n_levels);

        for(std::size_t i = 1; i < n_p1_levels; ++i) {
            auto m_i = meshes[i - 1]->clone();

            {
                libMesh::MeshRefinement mesh_refinement(*m_i);
                mesh_refinement.make_flags_parallel_consistent();
                mesh_refinement.uniformly_refine(1);
            }

            meshes[i] = std::move(m_i);
            use_interpolation_at_level_[i] = true;
        }

        if(is_second_order) {
            //p-multigrid does not need l2-projection but just interpolation
            use_interpolation_at_level_[n_p1_levels - 1] = use_interpolation_;
            use_interpolation_at_level_[n_p1_levels] = true;
        } else {
            use_interpolation_at_level_[n_coarse_spaces - 1] = use_interpolation_;
        }
    }

    std::unique_ptr<libMesh::MeshBase> SemiGeometricMultigrid::generate_box_mesh(const libMesh::MeshBase &fine_mesh, const std::size_t n_levels)
    {
        auto bb = bounding_box(fine_mesh);
        auto r = bb.max() - bb.min();

        const std::size_t n_coarse_spaces = n_levels - 1;
        const auto dim = fine_mesh.mesh_dimension();
        auto m = make_unique<libMesh::DistributedMesh>(fine_mesh.comm());

        switch(dim) {
            case 2:
            {
                const int n_segments = std::max(2, int( std::round( std::sqrt(fine_mesh.n_nodes() / std::pow(4, n_coarse_spaces)) ) ) );
                const double max_r = std::max(r(0), r(1));

                const double aspect_ratio_x = r(0)/max_r;
                const double aspect_ratio_y = r(1)/max_r;

                const int nx = std::max(int(std::round(n_segments * aspect_ratio_x)), 2);
                const int ny = std::max(int(std::round(n_segments * aspect_ratio_y)), 2);

                libMesh::MeshTools::Generation::build_square (
                    *m,
                    nx, ny,
                    bb.min()(0), bb.max()(0),
                    bb.min()(1), bb.max()(1),
                    libMesh::QUAD4);

                break;
            }

            case 3:
            {
                const int n_segments = std::max(2, int( std::round( std::cbrt(fine_mesh.n_nodes() / std::pow(8, n_coarse_spaces)) ) ) );

                const double max_r = std::max(r(0), std::max(r(1), r(2)));

                const double aspect_ratio_x = r(0)/max_r;
                const double aspect_ratio_y = r(1)/max_r;
                const double aspect_ratio_z = r(2)/max_r;

                const int nx = std::max(int(std::round(n_segments * aspect_ratio_x)), 2);
                const int ny = std::max(int(std::round(n_segments * aspect_ratio_y)), 2);
                const int nz = std::max(int(std::round(n_segments * aspect_ratio_z)), 2);

                libMesh::MeshTools::Generation::build_cube (
                    *m,
                    nx, ny, nz,
                    bb.min()(0), bb.max()(0),
                    bb.min()(1), bb.max()(1),
                    bb.min()(2), bb.max()(2),
                    libMesh::HEX8);

                break;
            }

            default:
            {
                assert(false && "implement me");
                break;
            }
        }

        return std::move(m);
    }

    void SemiGeometricMultigrid::init(const libMesh::System &es, const std::size_t n_levels)
    {
        const auto &mesh = es.get_mesh();
        const auto &dof_map = es.get_dof_map();
        const auto dim = mesh.mesh_dimension();

        bool is_second_order = dof_map.variable_type(0).order == 2;

        mg.fix_semidefinite_operators(true);

        generate_coarse_meshes(mesh, n_levels, dof_map.variable_type(0).order);

        const std::size_t n_coarse_spaces = n_levels - 1;
        equation_systems.resize(n_coarse_spaces);
        equation_systems[0] = std::make_shared<libMesh::EquationSystems>(*meshes[0]);

        for(std::size_t i = 0; i < n_coarse_spaces; ++i) {
            equation_systems[i] = std::make_shared<libMesh::EquationSystems>(*meshes[i]);

            auto &sys = equation_systems[i]->add_system<libMesh::LinearImplicitSystem>(es.name());

            for(unsigned int v = 0; v < dof_map.n_variables(); ++v) {
                const auto &var = dof_map.variable(v);
                //FIXME?
                sys.add_variable(var.name(), libMesh::FIRST, libMesh::LAGRANGE);
            }

            sys.init();

            std::cout << "level: " << i
                      << " n_dofs: " << equation_systems[i]->get_system(0).get_dof_map().n_dofs()
                      << " n_elems: " << equation_systems[i]->get_mesh().n_active_elem() << std::endl;
        }

        interpolators_.resize(n_coarse_spaces);
        moonolith::Communicator comm(mesh.comm().get());

        UVector d_diag;
        for(std::size_t i = 1; i < n_coarse_spaces; ++i) {
            interpolators_[i-1] = std::make_shared<USparseMatrix>();

            bool success = assemble_volume_transfer(
                comm,
                make_ref(*meshes[i-1]),
                make_ref(*meshes[i]),
                make_ref(equation_systems[i-1]->get_system(0).get_dof_map()),
                make_ref(equation_systems[i]->get_system(0).get_dof_map()),
                0,
                0,
                true,
                dof_map.n_variables(),
                *interpolators_[i-1],
                {},
                use_interpolation_at_level_[i - 1]
            ); assert(success);

            UVector d_diag;

            make_d(*interpolators_[i-1], d_diag);
            *interpolators_[i-1] = diag(1./d_diag) * *interpolators_[i-1];
        }

        interpolators_[n_coarse_spaces-1] = std::make_shared<USparseMatrix>();
        bool success = assemble_volume_transfer(
            comm,
            make_ref(*meshes[n_coarse_spaces-1]),
            make_ref(const_cast<libMesh::MeshBase &>(mesh)),
            make_ref(equation_systems[n_coarse_spaces-1]->get_system(0).get_dof_map()),
            make_ref(const_cast<libMesh::DofMap &>(dof_map)),
            0,
            0,
            true,
            dof_map.n_variables(),
            *interpolators_[n_coarse_spaces-1],
            {},
            use_interpolation_at_level_[n_coarse_spaces - 1]
        ); assert(success);

        make_d(*interpolators_[n_coarse_spaces-1], d_diag);
        *interpolators_[n_coarse_spaces-1] = diag(1./d_diag) * *interpolators_[n_coarse_spaces-1];

        if(mg.verbose()) {
            for(const auto &e : equation_systems) {
                std::cout << "dofs: " << e->get_system(0).get_dof_map().n_dofs() << std::endl;
            }

            std::cout << "dofs: " << es.get_dof_map().n_dofs() << std::endl;
        }

        mg.set_transfer_operators(interpolators_);
    }

    void SemiGeometricMultigrid::update(const std::shared_ptr<const USparseMatrix> &op)
    {
        mg.update(op);

#ifndef WITH_TRILINOS_ALGEBRA

        //hacky
        if(is_block_solver_) {
            for(SizeType i = 0; i < mg.n_levels(); ++i) {
                const_cast<USparseMatrix &>(mg.level(i).A()).convert_to_mat_baij(meshes[0]->mesh_dimension());
            }
        }

#endif //WITH_TRILINOS_ALGEBRA
    }

    bool SemiGeometricMultigrid::apply(const UVector &rhs, UVector &sol)
    {
        return mg.apply(rhs, sol);
    }

    void SemiGeometricMultigrid::update_contact(Contact &contact)
    {
        // const auto last_interp = mg.n_levels() - 2;
        // auto c_I = std::make_shared<USparseMatrix>();
        // *c_I = transpose(contact.complete_transformation) * *interpolators_[last_interp];

        
        // mg.update_transfer(last_interp, std::make_shared<MatrixTransfer<USparseMatrix, UVector>>(c_I));
    }
}


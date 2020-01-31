

#include "utopia_ContactSolver.hpp"

#include "utopia_ConvertToDiegoMeshApp.hpp"
#include "utopia_Newmark.hpp"
#include "utopia_LibMeshBackend.hpp"

#include "utopia_UIFunctionSpace.hpp"
#include "utopia_UIForcingFunction.hpp"
#include "utopia_UIMesh.hpp"

#include "utopia_ui.hpp"
#include "utopia_libmesh_DiegoMesh.hpp"
#include "libmesh/mesh_refinement.h"
#include "libmesh/exodusII_io.h"
#include "utopia.hpp"

namespace utopia {


    void ConvertToDiegoMeshApp::run(Input &in)
    {
        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));

        bool is_time_series = false;
        in.get("is-time-series", is_time_series);

        std::string folder = "./";
        in.get("folder", folder);

        std::set<int> blocks;
        std::cout << "Blocks?  "  << std::endl;

        in.get("blocks", [&](Input &in) {
            in.get_all([&](Input &in) {
                int block_id = -1;
                in.get("id", block_id);
                blocks.insert(block_id);
            });
        });

        std::cout << " ... " << blocks.size() << std::endl;

        if(!is_time_series) {
            in.get("mesh", mesh);
            in.get("space", space);

            if(!blocks.empty()) {
                std::cout << "Exporting blocks" << std::endl;
                DiegoMeshWriter().write(folder, mesh.mesh(), blocks);
                return;
            }

            DiegoMeshWriter().write(folder, mesh.mesh());

        } else {

            bool apply_disp = false;
            in.get("apply-disp", apply_disp);

            libMesh::ExodusII_IO io(mesh.mesh());

            in.get("mesh", [&](Input &in) {
                std::string path;
                in.get("path", path);
                io.read(path.c_str());
            });

            mesh.mesh().prepare_for_use();

            in.get("space", space);

            bool export_field_only = false;
            in.get("export-field-only", export_field_only);


            DiegoMeshWriter::Selector select;
            if(!export_field_only) {
                if(!blocks.empty()) {
                    std::cout << "Exporting blocks" << std::endl;
                    DiegoMeshWriter().write(folder, mesh.mesh(), blocks);
                    DiegoMeshWriter().build_index(mesh.mesh(), blocks, select);
                } else {
                    DiegoMeshWriter().write(folder, mesh.mesh());
                }
            }

            UVector data = local_zeros(space.space()[0].dof_map().n_local_dofs());
            UVector aux;

            int n_t = io.get_num_time_steps();

            int from = 0;
            int to   = n_t;
            int every = 1;

            in.get("from", from);
            in.get("to", to);
            in.get("every", every);

            to = std::min(n_t, to);

            for(uint i = from; i < to; i += every) {
                space.space().each([&](const int subspace_id, LibMeshFunctionSpace &s) {
                    const auto &name = s.var_name();
                    // std::cout << "reading var \"" << name << "\"" << std::endl;
                    io.copy_nodal_solution(
                        s.equation_system(),
                        name,
                        name,
                        i+1
                    );
                });

                utopia::convert(*space.space()[0].equation_system().solution, data);

                if(apply_disp) {
                    UScalar norm_disp = norm2(data);
                    std::cout << "norm2(disp) = " << norm_disp << std::endl;
                    deform_mesh(
                            mesh.mesh(),
                            space.space()[0].dof_map(),
                            data
                    );

                    std::string step_name = "t" + std::to_string(i + 1) + "_";

                    if(!blocks.empty()) {
                        DiegoMeshWriter().write_coords(folder, mesh.mesh(), select, step_name);
                    } else {
                        DiegoMeshWriter().write_coords(folder, mesh.mesh(), step_name);
                    }

                    aux = -data;

                    deform_mesh(
                        mesh.mesh(),
                        space.space()[0].dof_map(),
                        aux
                    );
                } else {

                    space.space().each([&](const int subspace_id, LibMeshFunctionSpace &s) {
                        const auto &name = s.var_name();

                        DiegoMeshWriter().write_field(
                            folder,
                            name + "_t" + std::to_string(i + 1),
                            mesh.mesh(),
                            data,
                            subspace_id
                        );
                    });
                }

            }
        }
    }
}


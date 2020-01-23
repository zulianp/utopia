

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

namespace utopia {


    void ConvertToDiegoMeshApp::run(Input &in)
    {
        UIMesh<libMesh::DistributedMesh> mesh(comm());
        UIFunctionSpace<LibMeshFunctionSpace> space(make_ref(mesh));
        in.get("mesh", mesh);
        in.get("space", space);

        std::string folder = "./";
        in.get("folder", folder);

        DiegoMeshWriter().write(folder, mesh.mesh());
    }
}


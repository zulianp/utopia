// #include "utopia_Main.hpp"
// #include "utopia_stk_Mesh.hpp"

// int main(const int argc, char *argv[]) { return UTOPIA_MAIN(argc, argv); }

// using namespace utopia;

// void stk_mesh_app(Input &in) {
//     using Mesh = utopia::stk::Mesh;

//     Mesh mesh;
//     mesh.read("dummy");
// }

// UTOPIA_REGISTER_APP(stk_mesh_app);

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>

int main(const int argc, char *argv[]) {
    using IOBroker = ::stk::io::StkMeshIoBroker;
    IOBroker broker(MPI_COMM_WORLD);

    return 0;
}
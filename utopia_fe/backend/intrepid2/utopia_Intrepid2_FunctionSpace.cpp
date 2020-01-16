#include "utopia_Intrepid2_FunctionSpace.hpp"
//#include "libmesh/nemesis_io.h"
#include "utopia_Intrepid2Backend.hpp"

namespace utopia {

    void write(const Path &path, Intrepid2FunctionSpace &space, TUSerialVector &x)
    {
       // utopia::convert(x, *space.equation_system().solution);
       // space.equation_system().solution->close();
  //      libMesh::Nemesis_IO(space.mesh()).write_equation_systems(path.to_string(), space.equation_systems());
    }
}

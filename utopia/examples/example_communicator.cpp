#include "utopia.hpp"
#include "utopia_Version.hpp"

#include <iostream>
#include <sstream>

template <class Vector>
void example_communicator(utopia::Input &in) {
    using Comm = typename utopia::Traits<Vector>::Communicator;

    // Equivalent to MPI comm world (serial for serial backends)
    auto &&world = Comm::world();

    // Serial communicator
    auto &&self = Comm::self();

    // Typically same as world
    auto &&default_comm = Comm::get_default();

    std::stringstream ss;

    std::string msg = "Comm sizes:";

    // You can change it from the command-line by using the `-msg <Your message> option`
    in.get("msg", msg);

    ss << msg << "\n";
    ss << "World:   " << world.size() << "\n";
    ss << "Self:    " << self.size() << "\n";
    ss << "Default: " << default_comm.size() << "\n";
    ss << "Default == World: " << (default_comm.same(world) ? "Yes" : "No") << std::endl;

    // Synchronized print to the standard output
    world.synched_print(ss.str(), std::cout);
}

// Run it with `./examples/example_communicator`
int main(int argc, char **argv) {
    using namespace utopia;

    // #ifdef WITH_PETSC
    //     using VectorT = PetscVector;
    // #else
    // #ifdef WITH_TRILINOS
    using VectorT = TpetraVectord;
    // #else
    //     using VectorT = BlasVectord;
    // #endif
    // #endif

    Utopia::Init(argc, argv);

    InputParameters params;
    params.init(argc, argv);

    example_communicator<VectorT>(params);

    return Utopia::Finalize();
}
